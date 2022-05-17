sisglmlss <- function(x, y, family = gamlss.dist::NO(), tune = c("BIC", "HQC", "AIC", "CV"), nfolds = 10L,
                      k.ses = NULL, nsis = NULL, iter.max = 10L, varISIS = c("van", "aggr"),
                      q = 1, seed = NULL, standardize = TRUE, ISIS.trace = FALSE,
                      gamlss.o.control = gamlss::gamlss.control(trace = FALSE, save = FALSE),
                      gamlss.i.control = gamlss::glim.control(),
                      gnet.control = gamlss.lasso2::gnet.control(nlambda = 90L, standardize = standardize),
                      parallel = FALSE, ncores = 2L, ...){
  
  tcall <- match.call()
  tune <- match.arg(tune)
  varISIS <- match.arg(varISIS)
  
  if(!is.matrix(x)) stop("x has to be a numeric matrix.")
  if(nrow(x) != length(y)) stop("x and y are not of compatible dimension.")
  if(!is.numeric(nfolds)) stop("nfolds has to be numeric.")
  if(!inherits(family, "gamlss.family")) stop("family has to be a gamlss.family object.")
  if(!is.null(nsis) && !is.numeric(nsis)) stop("nsis has to be numeric.")
  
  fit <- sisglmlss.fit(x, y, family, tune, nfolds, k.ses, nsis, iter.max, varISIS, q, seed, standardize, ISIS.trace,
                       gamlss.o.control, gamlss.i.control, gnet.control, parallel, ncores)
  fit$call <- tcall
  class(fit) <- "sisglmlss"
  fit
}



################################# FITTING FUNCTIONS #################################

sisglmlss.fit <- function(x, y, family, tune, nfolds, k.ses, nsis, iter.max, varISIS, q, seed, standardize, ISIS.trace,
                          gamlss.o.control, gamlss.i.control, gnet.control, parallel, ncores){
  n <- nrow(x)
  p <- ncol(x)
  if(is.null(k.ses)){
    k.ses <- c(1, rep(0, family$nopar - 1))
  } else {
    stopifnot(length(k.ses) == family$nopar)
  }
  if(!is.null(nsis)){
    if(length(nsis) != family$nopar) stop("For local selection, nsis needs to be a vector of length ", family$nopar, ".")
  } else {
    nsis <- rep(floor(compute.nsis(family, "van", n, p) / family$nopar), family$nopar)
  }
  names(nsis) <- names(family$parameters)
  
  if(isTRUE(standardize)){
    x <- scale(x)
  }
  if(!is.null(seed)){
    set.seed(seed)
  }
  #Create bootstrap-based sample splits if varISIS is not "van":
  if(varISIS == "van"){
    bootsplits <- NULL
  } else {
    bootsplits <- list(sample.int(n, n, replace = TRUE))
    alreadysampled <- unique(bootsplits[[1]])
    bootsplits[[2]] <- sample(c(setdiff(1:n, alreadysampled), sample.int(n, length(alreadysampled), replace = TRUE)))
  }
  # For local screening, we always need a data.frame containing the current selection (for code compactness and efficiency reasons):
  init.scr <- data.frame(var.ind = numeric(0), par = character(0))
  sels.allits <- list() #List to save all models for later comparison
  iter.id <- 1L
  for(par in names(family$parameters)){
    cat("### PARAMETER ", toupper(par), " ###\n")
    if(isTRUE(ISIS.trace)) cat("### Iteration ", iter.id, "\n")
    init.scr <- init.par.scr(x = x, y = y, bootsplits = bootsplits, par = par, sel.other.pars = init.scr, family = family, nsis = nsis, varISIS = varISIS,
                             q = q, gamlss.o.control = gamlss.o.control, gamlss.i.control = gamlss.i.control, parallel = parallel, ncores = ncores)
    init.scr <- check.singlevar(init.scr, p = p) # Single-variable check (can happen here for parameters other than mu)
    pen.sel <- pen.step(x = x, y = y, sel.mat = init.scr, family = family, tune = tune, nfolds = nfolds, k.ses = k.ses,
                        gnet.control = gnet.control)
    scr.sel <- pen.sel$sel
    if(isTRUE(ISIS.trace)) print(scr.sel)
    # Fully specify an equation for this parameter first, i.e. iterate only for par:
    repeat{
      sels.allits[[iter.id]] <- pen.sel$sel
      if(iter.id > 1L){
        # Check for first break condition: has the current selection already been found in any previous iteration, i.e.
        m.ids <- lapply(sels.allits[1:(iter.id - 1L)], identify.sel)
        any.equiv <- any(sapply(m.ids, setequal, y = identify.sel(pen.sel$sel)))
        if(any.equiv){
          cat("Model already selected\n")
          break
        }
      }
      iter.id <- iter.id + 1L
      if(isTRUE(ISIS.trace)) cat("### Iteration ", iter.id, "\n")
      scr.sel <- cond.par.scr(x = x, y = y, bootsplits = bootsplits, par = par, cur.sel = scr.sel, family = family, nsis = nsis, varISIS = varISIS, q = q,
                              gamlss.o.control = gamlss.o.control, gamlss.i.control = gamlss.i.control, parallel = parallel, ncores = ncores)
      # Check for second break condition: has the conditional screening step delivered any new variables?
      if(setequal(identify.sel(pen.sel$sel), identify.sel(scr.sel))){
        cat("Current screening step did not give any new variables.\n")
        break
      }
      scr.sel <- check.singlevar(scr.sel, p = p) # Single-variable check
      pen.sel <- pen.step(x = x, y = y, sel.mat = scr.sel, family = family, tune = tune, nfolds = nfolds, k.ses = k.ses,
                          gnet.control = gnet.control)
      if(isTRUE(ISIS.trace)) print(pen.sel$sel)
      # Check for next three break conditions (no variable survived pen.step, at least nsis variables selected or iter.max reached):
      if(nrow(pen.sel$sel) == 0L){
        cat("No variables remaining after penalized likelihood step in iteration", iter.id, ".\n")
        cat("Run again to try a different sample split (if varISIS %in% c('aggr', 'cons')) or use a more conservative variable screening approach!\n")
        break
      }
      if(all(table(factor(pen.sel$sel$par, levels = names(nsis))) >= nsis)){
        cat("Maximum number of variables (", nsis, ") selected for all distribution parameters.\n", sep = "")
        break
      }
      if(iter.id >= iter.max){
        cat("Maximum number of iterations (", iter.max, ") reached.\n", sep = "")
        break
      }
      scr.sel <- pen.sel$sel
    }
    init.scr <- pen.sel$sel
  }
  
  # Coefficient, intercept, lambda and selection extraction from (final (if ISIS)) penalized likelihood step:
  coef.beta <- pen.sel$betas
  a0 <- pen.sel$intercepts
  opt.lambdas <- pen.sel$opt.lambdas
  fin.sel <- pen.sel$sel
  for(i in seq_along(coef.beta)){
    coef.beta[[i]] <- c(a0[[i]], coef.beta[[i]])
    names(coef.beta[[i]])[1L] <- "(Intercept)"
  }
  return(list(init.scr = init.scr, fin.sel = fin.sel, coef.est = coef.beta, opt.lambdas = opt.lambdas,
              tuned.by = tune, fam = family))
}


init.par.scr <- function(x, y, bootsplits, par = c("mu", "sigma", "nu", "tau"), sel.other.pars, family, nsis, varISIS, q,
                         gamlss.o.control, gamlss.i.control, parallel, ncores){
  par <- match.arg(par)
  stopifnot(par %in% names(nsis))
  # ISIS implementation, heuristically compute k1 as in the global case:
  k1 <- floor(2*nsis / 3)
  if(varISIS == "van"){
    if(!is.null(bootsplits)) warning("Sample splits in bootsplits will be ignored if varISIS == 'van'.")
    par.devs <- compute.deviances(x, y, pars = par, cur.sel = sel.other.pars, family = family, gamlss.o.control = gamlss.o.control,
                                  gamlss.i.control = gamlss.i.control, parallel = parallel, ncores = ncores)
    dev.ranks <- rank(par.devs)
    counter <- 0L
    # Give it five different tries (permutations) to find relevant variables:
    repeat{
      counter <- counter + 1L
      rand.devs <- compute.deviances(x, y, pars = par, cur.sel = sel.other.pars, null.mod = TRUE, family = family,
                                     gamlss.o.control = gamlss.o.control, gamlss.i.control = gamlss.i.control,
                                     parallel = parallel, ncores = ncores)
      selection <- which(par.devs <= quantile(rand.devs, 1 - q))
      if(length(selection) > 0L) break
      if(counter >= 5L) break
    }
    perm.p <- length(selection)
    # Select at least 2, but at most k1 variables above the quantile threshold (and the highest below if above don't exist):
    selection <- which(dev.ranks <= max(min(perm.p, k1[par]), 2L))
  } else {
    if(is.null(bootsplits)) stop("Cannot perform varISIS 'aggr' or 'cons' without sample splits!")
    par.devs.s1 <- compute.deviances(x[bootsplits[[1L]], ], y[bootsplits[[1L]]], pars = par, cur.sel = sel.other.pars, family = family,
                                     gamlss.o.control = gamlss.o.control, gamlss.i.control = gamlss.i.control, parallel = parallel, ncores = ncores)
    par.devs.s2 <- compute.deviances(x[bootsplits[[2L]], ], y[bootsplits[[2L]]], pars = par, cur.sel = sel.other.pars, family = family,
                                     gamlss.o.control = gamlss.o.control, gamlss.i.control = gamlss.i.control, parallel = parallel, ncores = ncores)
    dev.ranks.s1 <- rank(par.devs.s1)
    dev.ranks.s2 <- rank(par.devs.s2)
    counter <- 0L
    # Give it five different tries (permutations) to find relevant variables:
    repeat{
      counter <- counter + 1L
      rand.devs.s1 <- compute.deviances(x[bootsplits[[1L]], ], y[bootsplits[[1L]]], pars = par, cur.sel = sel.other.pars, null.mod = TRUE, family = family,
                                        gamlss.o.control = gamlss.o.control, gamlss.i.control = gamlss.i.control, parallel = parallel, ncores = ncores)
      rand.devs.s2 <- compute.deviances(x[bootsplits[[2L]], ], y[bootsplits[[2L]]], pars = par, cur.sel = sel.other.pars, null.mod = TRUE, family = family,
                                        gamlss.o.control = gamlss.o.control, gamlss.i.control = gamlss.i.control, parallel = parallel, ncores = ncores)
      selection.s1 <- which(par.devs.s1 <= quantile(rand.devs.s1, 1 - q))
      selection.s2 <- which(par.devs.s2 <= quantile(rand.devs.s2, 1 - q))
      # Break only if we find variables breaking the deviance quantile threshold in both splits:
      if(length(selection.s1) > 0L && length(selection.s2) > 0L) break
      if(counter >= 5L) break
    }
    perm.p.s1 <- length(selection.s1)
    perm.p.s2 <- length(selection.s2)
    selection <- intersect(which(dev.ranks.s1 <= min(perm.p.s1, k1[par])),
                           which(dev.ranks.s2 <= min(perm.p.s2, k1[par])))
    if(length(selection) <= 1L){
      m <- compute.m(dev.ranks.s1, dev.ranks.s2, common = 2L, min = (k1[par] + 1L))
      selection <- intersect(which(dev.ranks.s1 <= m), which(dev.ranks.s2 <= m))
    }
  }
  selection <- data.frame(var.ind = selection, par = par)
  rbind(sel.other.pars, selection)
}


cond.par.scr <- function(x, y, bootsplits, par = c("mu", "sigma", "nu", "tau"), cur.sel, family, nsis, varISIS, q,
                         gamlss.o.control, gamlss.i.control, parallel, ncores){
  par <- match.arg(par)
  stopifnot(par %in% names(nsis))
  p <- ncol(x)
  kl <- nsis - table(factor(cur.sel$par, levels = names(nsis))) #Remaining variables (parameter-specific)
  if(kl[par] == 0L) return(cur.sel)
  if(varISIS == "van"){
    if(!is.null(bootsplits)) warning("Sample splits in bootsplits will be ignored if varISIS == 'van'.")
    par.devs <- compute.deviances(x, y, pars = par, cur.sel = cur.sel, family = family, gamlss.o.control = gamlss.o.control,
                                  gamlss.i.control = gamlss.i.control, parallel = parallel, ncores = ncores)
    dev.ranks <- rank(par.devs)
    rand.devs <- compute.deviances(x, y, pars = par, cur.sel = cur.sel, null.mod = TRUE, family = family,
                                   gamlss.o.control = gamlss.o.control, gamlss.i.control = gamlss.i.control,
                                   parallel = parallel, ncores = ncores)
    new.sel <- which(par.devs <= quantile(rand.devs, 1 - q, na.rm = TRUE))
    new.sel <- which(dev.ranks <= min(length(new.sel), kl[par]))
  } else {
    if(is.null(bootsplits)) stop("Cannot perform varISIS 'aggr' or 'cons' without sample splits!")
    par.devs.s1 <- compute.deviances(x[bootsplits[[1L]], ], y[bootsplits[[1L]]], pars = par, cur.sel = cur.sel, family = family,
                                     gamlss.o.control = gamlss.o.control, gamlss.i.control = gamlss.i.control, parallel = parallel, ncores = ncores)
    par.devs.s2 <- compute.deviances(x[bootsplits[[2L]], ], y[bootsplits[[2L]]], pars = par, cur.sel = cur.sel, family = family,
                                     gamlss.o.control = gamlss.o.control, gamlss.i.control = gamlss.i.control, parallel = parallel, ncores = ncores)
    dev.ranks.s1 <- rank(par.devs.s1)
    dev.ranks.s2 <- rank(par.devs.s2)
    #Find distribution of conditional marginal deviances under a null model in each bootstrap sample:
    rand.devs.s1 <- compute.deviances(x[bootsplits[[1L]], ], y[bootsplits[[1L]]], pars = par, cur.sel = cur.sel, null.mod = TRUE, family = family,
                                      gamlss.o.control = gamlss.o.control, gamlss.i.control = gamlss.i.control, parallel = parallel, ncores = ncores)
    rand.devs.s2 <- compute.deviances(x[bootsplits[[2L]], ], y[bootsplits[[2L]]], pars = par, cur.sel = cur.sel, null.mod = TRUE, family = family,
                                      gamlss.o.control = gamlss.o.control, gamlss.i.control = gamlss.i.control, parallel = parallel, ncores = ncores)
    selection.s1 <- which(par.devs.s1 <= quantile(rand.devs.s1, 1 - q, na.rm = TRUE))
    selection.s2 <- which(par.devs.s2 <= quantile(rand.devs.s2, 1 - q, na.rm = TRUE))
    perm.p.s1 <- length(selection.s1)
    perm.p.s2 <- length(selection.s2)
    new.sel <- intersect(which(dev.ranks.s1 <= min(perm.p.s1, kl[par])),
                         which(dev.ranks.s2 <= min(perm.p.s2, kl[par])))
  }
  new.sel <- data.frame(var.ind = new.sel, par = rep(par, length(new.sel))) #can be of nrow = 0, if no new vars were found!
  selection <- rbind(cur.sel, new.sel)
  selection <- selection[order(selection$par, selection$var.ind), ]
  rownames(selection) <- as.character(seq_len(nrow(selection)))
  selection
}


################################# HELPER FUNCTIONS RELEVANT FOR BOTH SEL.TYPES #################################

identify.sel <- function(sel.mat){
  paste(sel.mat$par, sel.mat$var.ind, sep = "_")
}

current.distlist <- function(){
  dists <- c("BB", "BCCG", "BCCGo", "BCCGuntr", "BCPE", "BCPEo", "BCPEuntr", "BCT", "BCTo", "BCTuntr", "BE", "BEINF", "BEINF0",
             "BEINF1", "BEo", "BEOI", "BEZI", "BI", "BNB", "DBI", "DBURR12", "DEL", "DPO", "DPO1", "EGB2", "EXP", "GA", "GAF",
             "GB1", "GB2", "GEOM", "GEOMo", "GG", "GIG", "GP", "GPO", "GT", "GU", "IG", "IGAMMA", "JSU", "JSUo", "LG", "LNO",
             "LO", "LOGITNO", "LOGNO", "LOGNO2", "LQNO", "MN3", "MN4", "MN5", "MULTIN", "NBF", "NBI", "NBII", "NET", "NO",
             "NO2", "NOF", "PARETO", "PARETO1", "PARETO1o", "PARETO2", "PARETO2o", "PE", "PE2", "PIG", "PIG2", "PO", "RG",
             "RGE", "SEP", "SEP1", "SEP2", "SEP3", "SEP4", "SHASH", "SHASHo", "SHASHo2", "SI", "SICHEL", "SIMPLEX", "SN1",
             "SN2", "SST", "ST1", "ST2", "ST3", "ST3C", "ST4", "ST5", "TF", "TF2", "WARING", "WEI", "WEI2", "WEI3", "YULE",
             "ZABB", "ZABI", "ZABNB", "ZAGA", "ZAIG", "ZALG", "ZANBI", "ZAP", "ZAPIG", "ZASICHEL", "ZAZIPF", "ZIBB", "ZIBI",
             "ZIBNB", "ZINBF", "ZINBI", "ZIP", "ZIP2", "ZIPF", "ZIPIG", "ZISICHEL")
  # Categorization according to Fan et al. (2009) and Saldana and Feng (2018): model-based choises of nsis,
  # favouring smaller values in models where the response provides less information:
  # Category 1: unbounded real-valued response -> n/log(n)
  # Category 2: bounded real-valued response or unbounded integer-valued -> n/(2*log(n))
  # Category 3: bounded integer-valued -> n/(4*log(n))
  # Assign distributions to categories:
  cats <- c(3, 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 2, 2,
            2, 2, 2, 2, 3, 3, 2, 3, 2, 2, 2, 1, 2, 2, 2,
            2, 2, 3, 3, 2, 2, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2,
            1, 2, 2, 2, 1, 3, 3, 3, 3, 3, 3, 3, 1, 1,
            1, 1, 2, 2, 2, 2, 2, 1, 1, 2, 2, 2, 1,
            2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,
            3, 3, 3, 2, 2, 2, 3, 2, 2, 2, 2, 3, 3,
            3, 3, 3, 2, 2, 2, 2, 2)
  data.frame(dist = dists, cat = cats)
}

#TBD: Maybe adapt here, since we are selecting variables over K equations (i.e. select MORE than in the single-equation context)?
compute.nsis <- function(family, varISIS, n, p){
  if(varISIS == "aggr") return(floor(n/log(n)))
  if(p < n) return(p)
  family.name <- family$family[1L]
  dist.list <- current.distlist()
  dist.id <- match(family.name, dist.list$dist)
  if(is.na(dist.id)) stop("Family '", family.name, "' not found, please adapt distlist.")
  cat <- dist.list$cat[dist.id]
  if(cat == 1) return(floor(n/log(n)))
  if(cat == 2) return(floor(n/(2 * log(n))))
  return(floor(n/(4 * log(n))))
}

#Equivalent of SIS:::intensure
compute.m <- function(rank.mat1, rank.mat2, common = 2L, min = common){
  # Search to rank matrices, where they have at least "common" variables in common
  # (also depending on the regression equation of course). Start search at "min"
  # (which by default is "common") since we obviously need at least to look at the
  # top "common" ranks to find "common" variables in common:
  for(i in min:length(rank.mat1)){
    inters <- intersect(which(rank.mat1 <= i), which(rank.mat2 <= i))
    if(length(inters) >= common) return(i)
  }
}

compute.deviances <- function(x, y, pars = c("all", "mu", "sigma", "nu", "tau"), cur.sel = NULL, null.mod = FALSE, family,
                              gamlss.o.control, gamlss.i.control, parallel, ncores){
  n <- nrow(x)
  p <- ncol(x)
  # Local screening: compute marginal deviances by fitting univariate models for one regression parameter only.
  # NOTE: for this case, ALWAYS a cur.sel parameter has to be passed to the function. If the selection is for the
  # first parameter (mu) in the first iteration, then the data.frame should be of dimension c(0,2) with the right names
  # (those are var.ind and par).
  if(is.null(cur.sel)) stop("For local selection, a cur.sel argument is always required. If the model is empty so far, pass data.frame with 0 rows.")
  out <- matrix(NA, p, 1L)
  if(isTRUE(null.mod)){
    perm.int <- sample.int(n, n)
    basestring <- "x[perm.int, %d]"
  } else {
    basestring <- "x[, %d]"
  }
  cur.spec <- constr.formula.nopen(cur.sel, with.data = FALSE)
  cur.spec <- lapply(cur.spec, as.formula, env = environment())
  args_list <- c(cur.spec, list(family = family, control = gamlss.o.control, i.control = gamlss.i.control))
  form_name <- switch(pars, mu = "formula", paste0(pars, ".formula"))
  par_inds <- setdiff(1:p, cur.sel$var.ind[cur.sel$par == pars])
  if(isTRUE(parallel)){
    get.dev.par <- function(i, envir){
      args_list[[form_name]] <- as.formula(paste0("y ~ ", paste(c(sprintf("x[, %d]", cur.sel$var.ind[cur.sel$par == pars]),
                                                                  sprintf(basestring, i)), collapse = " + ")),
                                           env = envir)
      m <- try(do.call(gamlss::gamlss, args_list), silent = TRUE)
      #If RS() didn't work, try CG(), otherwise NA:
      if(inherits(m, "try-error")){
        args_listCG <- c(args_list, list(method = quote(CG())))
        m <- try(do.call(gamlss::gamlss, args_listCG), silent = TRUE)
      }
      if(inherits(m, "try-error")) return(NA_real_) else return(deviance(m))
    }
    eval.env <- environment()
    out[par_inds,1] <- unlist(parallel::mclapply(parallel:::splitList(par_inds, ncores),
                                                 function(p_inds) vapply(p_inds, get.dev.par, numeric(1), envir = eval.env),
                                                 mc.cores = ncores))
  } else {
    for(i in par_inds){
      args_list[[form_name]] <- as.formula(paste0("y ~ ", paste(c(sprintf("x[, %d]", cur.sel$var.ind[cur.sel$par == pars]),
                                                                  sprintf(basestring, i)), collapse = " + ")))
      m <- try(do.call(gamlss::gamlss, args_list), silent = TRUE)
      #If RS() didn't work, try CG(), otherwise NA:
      if(inherits(m, "try-error")){
        args_listCG <- c(args_list, list(method = quote(CG())))
        m <- try(do.call(gamlss::gamlss, args_listCG), silent = TRUE)
      }
      out[i,1] <- if(inherits(m, "try-error")) NA else deviance(m)
    }
  }
  out
}


#Equivalent of SIS::tune.fit function:
pen.step <- function(x, y, sel.mat, family, tune, nfolds, k.ses, gnet.control){
  n <- nrow(x)
  p <- ncol(x)
  d <- setNames(as.data.frame(cbind(x, y)), c(sprintf("V%d", 1:p), "y"))
  if(tune == "CV"){
    #Model tuning based on cross validation:
    CVfolds <- parallel:::splitList(sample.int(n, n), nfolds)
    CVfolds <- lapply(1:nfolds, function(i) do.call(c,CVfolds[(1:nfolds)[-i]]))
    pen.model.form <- constr.formula.pen(sel.mat, type = "CV", k.ses)
  } else {
    #Model tuning based on information criterion (set in "tune"):
    pen.model.form <- constr.formula.pen(sel.mat, type = "IC")
  }
  # We have to pass on the current function evaluation environment to as.formula as otherwise as.formula would bind
  # the temporary lapply-eval environment to the created formula objects, in which we would then not find the variables
  # of interest (i.e. the ones included in the gnet terms (such as sel.mat)) the ones for the model estimation:
  pen.model.form <- lapply(pen.model.form, as.formula, env = environment())
  #Important time and memory saver here: pass in only relevant columns of d:
  args_list <- c(pen.model.form, list(data = d[, sort(unique(sel.mat$var.ind))], family = family, control = gamlss::gamlss.control(trace = FALSE),
                                      i.control = gamlss::glim.control(cyc = 50L)))
  repeat{
    if(args_list$i.control$cyc < 1L){
      args_list$i.control$cyc <- 1L
      args_list$i.control$bf.cyc <- args_list$i.control$bf.cyc - 1L
      pen.fit <- try(do.call(gamlss::gamlss, args_list), silent = TRUE)
      if(!inherits(pen.fit, "try-error")){
        warning("Minimum number of i.control$cyc reached.")
        break
      }
      if(args_list$i.control$bf.cyc == 1L) stop("Minimum number of i.control$cyc and i.control$bf.cyc reached.")
    } else {
      pen.fit <- try(do.call(gamlss::gamlss, args_list), silent = TRUE)
      if(!inherits(pen.fit, "try-error")) break
      args_list$i.control$cyc <- args_list$i.control$cyc - 1L
    }
  }
  # The smooth terms that come back from estimating the gamlss with lasso are contained in a list that is accessible via
  # getSmo(pen.fit, "mu"), for example. The structure of this list differs, depending on whether tuning is based on IC or CV:
  # IC: 
  # With IC, the list of length 2. The first contains the object of class c("elnet", "glmnet") that represents the glmnet model
  # fit (lasso path) for a given parameter. The second contains the information criterion selected (in slot "OPTCRIT") and
  # the optimal chosen lambda (optlambda) that minimizes the IC, among other things. It also contains the regression coefficients
  # of the model contained in the first list element, specific to the optimal lambda.
  # CV:
  # With CV, the list is of length (nfolds + 2). The first nfolds elements of this list contain objects of class c("elnet", "glmnet")
  # that represent the glmnet model fits in each fold for a given parameter (e.g. "mu"). To explain these slots,
  # first consider what happens behind the scenes when calling gamlss::gamlss with a gnet regularization tuned with CV:
  # The lasso paths in the respective folds have all been fitted with the SAME lambda sequence thanks to an adaption in the
  # gamlss.lasso package implemented by me (wasn't the case before, see my comments in gamlss.lasso:::gamlss.gnet source
  # code). The out-of-fold (OOF) MSEs are then computed (and saved to be accessed via getSmo(pen.fit, "mu")[[nfolds + 2L]]$OPTCRIT)
  # and the best lambda is then the one that is the biggest among those, for which the avg. OOF MSE is at most as big as
  # the minimal OOF MSE (across all lambda) PLUS "k.se = 1" standard error(s). With the lambda chosen in this way,
  # the gamlss.lasso with type = "sel" then fits one more model to ALL data with the same lambda sequence as before. The
  # corresponding c("elnet", "glmnet") object is returned as the [[nfolds + 1L]] element of the list. Finally, the last
  # element of the list contains a summary of all the fits, namely, again the complete lambda sequence, all the OOF MSEs (see above),
  # the position of the best lambda in the sequence (optid), the regression coefficients (BETA) from the model estimated on the
  # whole data with the optimal lambda (so it is getSmo(pen.fit, "mu")[[nfolds + 1L]]$beta[,optid]) and the optimal lambda
  # (among other things).
  
  # We extract the optimal lambdas, betas and intercepts in each reg. equation and adapt our current selection matrix
  # based on which additional variables were kicked out by the lasso:
  pen.lambdas <- setNames(lapply(names(family$parameters), function(n) tail(getSmo(pen.fit,n), 1)[[1]]$optlambda),
                          names(family$parameters))
  pen.betas <- setNames(lapply(names(family$parameters), function(n) tail(getSmo(pen.fit,n), 1)[[1]]$beta),
                        names(family$parameters))
  pen.intercepts <- setNames(lapply(names(family$parameters), function(n) unname(coef(pen.fit, n)[1])),
                             names(family$parameters))
  pen.sel <- sel.mat[do.call(c, pen.betas) != 0, ]
  pen.betas <- lapply(pen.betas, function(x) x[x != 0])
  list(sel = pen.sel, intercepts = pen.intercepts, betas = pen.betas, opt.lambdas = pen.lambdas)
}

#Create model formulas for gamlss estimation with penalized likelihood
constr.formula.pen <- function(sel, type = c("IC","CV"), k.ses = NULL){
  if(!is.null(k.ses) && type == "IC") warning("k.ses only used with CV, will be ignored.")
  type <- match.arg(type)
  choices <- c("mu", "sigma", "nu", "tau")
  if(type == "CV"){
    form.out <- lapply(1:4, function(i){
      #1-SE rule by default (k.se = 1):
      if(any(sel$par == choices[i])) return(paste0("y ~ gnet(x.vars = sprintf(\"V%d\",sel.mat$var.ind[sel.mat$par == \"",
                                                   choices[i], "\"]), method = \"CV\", type = \"sel\", k.se = ",
                                                   k.ses[i], ", adaptive = NULL, subsets = CVfolds, control = gnet.control)"))
      return("y ~ 1")
    })
  } else {
    form.out <- lapply(1:4, function(i){
      if(any(sel$par == choices[i])) return(paste0("y ~ gnet(x.vars = sprintf(\"V%d\",sel.mat$var.ind[sel.mat$par == \"",
                                                   choices[i], "\"]), method = \"IC\", ICpen = tune, adaptive = NULL,
                                                   control = gnet.control)"))
      return("y ~ 1")
    })
  }
  names(form.out) <- c("formula", "sigma.formula", "nu.formula", "tau.formula")
  form.out
}

#Create model formulas for gamlss estimation with non-penalized likelihood:
constr.formula.nopen <- function(sel, with.data = TRUE){
  choices <- c("mu", "sigma", "nu", "tau")
  if(isTRUE(with.data)){
    form.out <- lapply(1:4, function(i){
      if(any(sel$par == choices[i])) return(paste0("y ~ ", paste(sprintf("V%d", sel$var.ind[sel$par == choices[i]]), collapse = " + ")))
      return("y ~ 1")
    })
  } else {
    form.out <- lapply(1:4, function(i){
      if(any(sel$par == choices[i])) return(paste0("y ~ ", paste(sprintf("x[, %d]", sel$var.ind[sel$par == choices[i]]), collapse = " + ")))
      return("y ~ 1")
    })
  }
  names(form.out) <- c("formula", "sigma.formula", "nu.formula", "tau.formula")
  form.out
}

#Check whether the current selection contains any distribution parameter, for which there is only one variable selected,
#which would not allow a penalized likelihood estimation:
check.singlevar <- function(cur.sel, p){
  par.table <- table(cur.sel$par)
  if(any(par.table == 1L)){
    add.pars <- lapply(names(par.table)[which(par.table == 1L)], function(n){
      id <- cur.sel$var.ind[cur.sel$par == n]
      data.frame(var.ind = p + 1 - id, par = n)
    })
    cur.sel <- do.call(rbind, c(list(cur.sel), add.pars))
    cur.sel <- cur.sel[order(cur.sel$par, cur.sel$var.ind), ]
  }
  cur.sel
}


##################### METHODS #####################
print.sisglmlss <- function(x, ...){
  cat("Family: ", x$fam$family[-1], "\n")
  cat("Call: ", deparse(x$call, width.cutoff = 50), "\n", fill = TRUE)
  cat("Final selection:\n")
  for(par in names(x$fam$parameters)){
    cat(par, ": ", paste(x$fin.sel$var.ind[x$fin.sel$par == par], collapse = ", "), "\n", sep = "")
  }
  cat("\n")
  for(par in names(x$fam$parameters)){
    cat(par, "coefficients:\n")
    print(x$coef.est[[par]])
  }
  invisible(x)
}
