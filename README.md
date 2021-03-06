
# sisglmlss

## Overview

**`sisglmlss`** is an R package to perform model selection via
(iterative) sure independence screening (ISIS) in generalized linear
models for location, scale and shape (GLMLSS). It implements the
algorithms `pvanISIS` and `pvarISIS` proposed in the paper by Amon and
Hornik (2022). For details, please see the original paper and the
package documentation.

## Installation

The source code is hosted on
[github](https://github.com/jamon-R/sisglmlss). To install the *R*
package, run the following code:

``` r
devtools::install_github("jamon-R/sisglmlss")
```

The package requires functionality from the packages `gamlss`,
`gamlss.dist` and `gamlss.lasso2`. While the former two are available on
[CRAN](https://cran.r-project.org/), the latter package is only
available in a version called `gamlss.lasso`, which does not perform
tuning based on cross-validation in a correct manner. We therefore had
to adapt the source code from the corresponding package ourselves. The
corresponding adaptation is also available on
[github](https://github.com/jamon-R/gamlss.lasso2) and has to be
installed prior to the installation of `sisglmlss`.

## Usage

The main work horse of the package is the function `sisglmlss`, which
allows to perform different variants of sure independence screening in
GLMLSSs. The function offers a lot of flexibility for the user.
Specifically, the user can decide, which family of distributions to use
for the response variable, how to tune the shrinkage parameter in the
LASSO step, whether to perform the `pvanISIS` or `pvarISIS` variant,
which quantile threshold to use for the data-driven likelihood
thresholding, whether or not to standardize the variables and how many
cores to use for the marginal maximum likelihood computations. Moreover,
the user also has full control of the inner workings of the model
fitting, as there are arguments to pass in output from the outer and
inner control functions of the
[gamlss](https://cran.r-project.org/web/packages/gamlss/index.html)
package.
