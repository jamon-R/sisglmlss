\name{sisglmlss}
\alias{sisglmlss}

\title{
Variable selection with iterative sure independence screening (ISIS) in GLMLSS
}
\description{
This function performs iterative sure independence screening (ISIS) in a generalized linear model for location, scale and shape. Based on predictor matrix \code{x} and response vector \code{y}, the function performs alternating marginal likelihood ratio screening and penalized likelihood steps for all distribution parameters of the response distribution specified via \code{family}. The function offers all kinds of detailed tuning options.}
\usage{
sisglmlss(x,
  y,
  family = gamlss.dist::NO(),
  tune = c("BIC", "HQC", "AIC", "CV"),
  nfolds = 10L,
  k.ses = NULL,
  nsis = NULL,
  iter.max = 10L,
  varISIS = c("van", "aggr"),
  q = 1,
  seed = NULL,
  standardize = TRUE,
  ISIS.trace = FALSE,
  gamlss.o.control = gamlss::gamlss.control(trace = FALSE, save = FALSE),
  gamlss.i.control = gamlss::glim.control(),
  gnet.control = gamlss.lasso2::gnet.control(nlambda = 90L, standardize = standardize),
  parallel = FALSE,
  ncores = 2L, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ The design matrix of dimension \code{n * p}, without an intercept. Each row is an observation vector.}
  \item{y}{ The response vector of dimension \code{n}. Has to be a valid response vector under the response distribution specified under \code{family}.}
  \item{family}{ A \code{\link[gamlss.dist]{gamlss.family}} object, which is used to define the distribution and the link functions of the various parameters. The supported distribution families can be found in \code{\link[gamlss.dist]{gamlss.family}}. }
  \item{tune}{ Method for tuning the regularization parameter of the penalized likelihood steps. Options are \code{"BIC"}, \code{"HQC"}, \code{"AIC"} and \code{"CV"} (k-fold cross validation). Default is \code{"BIC"}. }
  \item{nfolds}{ Number of folds to use in k-fold cross validation when \code{tune = "CV"}, the default is 10. }
  \item{k.ses}{  This parameter establishes how many times the standard deviation is summed to the mean to select the optimal lambda, when tuning the optimal regularization parameter via cross-validation, i.e. if \code{tune} is \code{"CV"}. By default, it is 1 for distribution parameter \code{"mu"} and 0 for all others. }
  \item{nsis}{ Maximum number of predictors to be used in each distribution parameter. Has to be either \code{NULL} or a vector of length \code{family$nopar}. If \code{NULL} (default), then the maximum number of equation-wise predictors is chosen in accordance to a matching table depending on the information contained in the response. See 'Details'. }
  \item{iter.max}{ Maximum number of ISIS iterations. Default is 10.}
  \item{varISIS}{ Specifies whether to perform vanilla-ISIS or aggr-ISIS. Vanilla-ISIS uses the entire sample provided to the function for screening and penalized likelihood steps. Aggr-ISIS creates two bootstrap samples of length (n) (i.e. draws randomly from the rows of \code{x} and \code{y} with replacement) and then keeps only variables are detected in the screening step in both bootstrap samples. It is more aggressive because it decreases the false discovery rate, but runs the risk of excluding signal. }
  \item{q}{ Quantile for calculating the data-driven threshold in the screening step. Marginal likelihoods are drawn from a null model by permuting the rows of \code{x} whilst keeping \code{y} unchanging, thereby decoupling response from predictors. From the null distribution of marginal likelihoods under this model, a quantile threshold is computed that the original variables have to break, in order to be considered for the model. \code{q} is thus between 0 and 1 and should be chosen close to 1. Default is 1. }
  \item{seed}{ Random seed to be set before generating any pseudo-random numbers. }
  \item{standardize}{ Logical indicating whether the columns of \code{x} should be standardized, i.e. subtract mean and divide by standard deviation. }
  \item{ISIS.trace}{ Logical indicating whether the current state of the variable selection process should be printed to the console at each iteration. Defaults to \code{FALSE}. }
  \item{gamlss.o.control}{ This sets the control parameters of the outer iterations algorithm using the function \code{\link[gamlss]{gamlss.control}}. }
  \item{gamlss.i.control}{ This sets the control parameters of the inner iterations of the RS algorithm using the function \code{\link[gamlss]{glim.control}}. }
  \item{gnet.control}{ This sets the control paramters for glmnet for the penalized likelihood step. To be set using the function \code{\link[gamlss.lasso2]{gnet.control}}. }
  \item{parallel}{ Specifies whether the marginal likelihood computations should be parallelized using parallel::mclapply. Strongly recommended for \code{p > 1000}, yet only possible on Mac and Linux systems. Defaults to \code{FALSE}. }
  \item{ncores}{ Specifies the number of cores to use for parallelization if \code{parallel} is \code{TRUE}. Should ideally be an integer multiple of the number of distribution parameters, otherwise there will be a warning. Defaults to 2. }
  \item{\dots}{ Extra arguments to be passed on to helper functions. Currently not in use.}
}
\details{
Put matching table from family to nsis here.
}
\value{
Returns an object of class \code{"sisglmlss"} with the following components: 
\item{init.scr }{Initial screening for the final parameter.}
\item{fin.sel }{Data.Frame containing the variable ids of all relevant variables for each distribution parameter.}
\item{coef.est }{Regression coefficient estimates from the final penalized likelihood step.}
\item{opt.lambdas }{Optimal regularization parameter for each distribution parameter in the final penalized likelihood step.}
\item{tuned.by }{Method of choosing the optimal regularization paramter, equal to input \code{tune}.}
\item{fam }{A \code{\link[gamlss.dist]{gamlss.family}} object specifying the employed response distribution.}
\item{call }{The call made to the function.}
}
\references{
Diego Franco Saldana and Yang Feng (2018) SIS: An R package for Sure Independence Screening in
Ultrahigh Dimensional Statistical Models, \emph{Journal of Statistical Software}, \bold{83}, 2, 1-25.

Jianqing Fan and Rui Song (2010) Sure Independence Screening in Generalized
Linear Models with NP-Dimensionality.  \emph{The Annals of Statistics},
\bold{38}, 3567-3604.

Jianqing Fan, Richard Samworth, and Yichao Wu (2009) Ultrahigh Dimensional
Feature Selection: Beyond the Linear Model. \emph{Journal of Machine
Learning Research}, \bold{10}, 2013-2038.

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{https://www.jstatsoft.org/v23/i07/}.
}
\author{ Julian Amon \email{julian.amon@wu.ac.at} }
\note{The following generic functions can be used with a sisglmlss object: \code{print}. }
\examples{
set.seed(123)
n <- 50L
p <- 100L
betas <- list(mu = c(1, 1.5, -0.5, -0.75, -1.25, rep(0, p - 5L)),
              sigma = c(0, 0, 0, -0.4, -0.2, 0.2, 0.3, 0.1, rep(0, p - 8L)))
x <- matrix(rnorm(n*p), n, p)
y <- rnorm(n, mean = colSums(t(x)*betas$mu), sd = exp(colSums(t(x)*betas$sigma)))
\dontrun{
    sel <- sisglmlss(x, y)
    print(sel)
}
}

\keyword{regression}
