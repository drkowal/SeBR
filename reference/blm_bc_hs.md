# Bayesian linear model with a Box-Cox transformation and a horseshoe prior

MCMC sampling for Bayesian linear regression with 1) a (known or
unknown) Box-Cox transformation and 2) a horseshoe prior for the
(possibly high-dimensional) regression coefficients.

## Usage

``` r
blm_bc_hs(
  y,
  X,
  X_test = X,
  lambda = NULL,
  sample_lambda = TRUE,
  only_theta = FALSE,
  nsave = 1000,
  nburn = 1000,
  nskip = 0,
  verbose = TRUE
)
```

## Arguments

- y:

  `n x 1` vector of observed counts

- X:

  `n x p` matrix of predictors (no intercept)

- X_test:

  `n_test x p` matrix of predictors for test data; default is the
  observed covariates `X`

- lambda:

  Box-Cox transformation; if NULL, estimate this parameter

- sample_lambda:

  logical; if TRUE, sample lambda, otherwise use the fixed value of
  lambda above or the MLE (if lambda unspecified)

- only_theta:

  logical; if TRUE, only return posterior draws of the regression
  coefficients (for speed)

- nsave:

  number of MCMC iterations to save

- nburn:

  number of MCMC iterations to discard

- nskip:

  number of MCMC iterations to skip between saving iterations, i.e.,
  save every (nskip + 1)th draw

- verbose:

  logical; if TRUE, print time remaining

## Value

a list with the following elements:

- `coefficients` the posterior mean of the regression coefficients

- `fitted.values` the posterior predictive mean at the test points
  `X_test`

- `post_theta`: `nsave x p` samples from the posterior distribution of
  the regression coefficients

- `post_ypred`: `nsave x n_test` samples from the posterior predictive
  distribution at test points `X_test`

- `post_g`: `nsave` posterior samples of the transformation evaluated at
  the unique `y` values

- `post_lambda`: `nsave` posterior samples of lambda

- `post_sigma`: `nsave` posterior samples of sigma

- `model`: the model fit (here, `blm_bc_hs`)

as well as the arguments passed in.

## Details

This function provides fully Bayesian inference for a transformed linear
model via MCMC sampling. The transformation is parametric from the
Box-Cox family, which has one parameter `lambda`. That parameter may be
fixed in advanced or learned from the data.

The horseshoe prior is especially useful for high-dimensional settings
with many (possibly correlated) covariates. This function uses a fast
Cholesky-forward/backward sampler when `p < n` and the Bhattacharya et
al. (\<https://doi.org/10.1093/biomet/asw042\>) sampler when `p > n`.
Thus, the sampler can scale linear in `n` (for fixed/small `p`) or
linear in `p` (for fixed/small `n`).

## Note

Box-Cox transformations may be useful in some cases, but in general we
recommend the nonparametric transformation in
[`sblm_hs`](https://drkowal.github.io/SeBR/reference/sblm_hs.md).

An intercept is automatically added to `X` and `X_test`. The
coefficients reported do \*not\* include this intercept parameter, since
it is not identified under more general transformation models (e.g.,
[`sblm_hs`](https://drkowal.github.io/SeBR/reference/sblm_hs.md)).

## Examples

``` r
# Simulate data from a transformed (sparse) linear model:
dat = simulate_tlm(n = 100, p = 50, g_type = 'step', prop_sig = 0.1)
y = dat$y; X = dat$X # training data

hist(y, breaks = 25) # marginal distribution


# Fit the Bayesian linear model with a Box-Cox transformation & a horseshoe prior:
fit = blm_bc_hs(y = y, X = X, verbose = FALSE)
names(fit) # what is returned
#>  [1] "coefficients"  "fitted.values" "post_theta"    "post_ypred"   
#>  [5] "post_g"        "post_lambda"   "post_sigma"    "model"        
#>  [9] "y"             "X"             "X_test"        "sample_lambda"
#> [13] "only_theta"   
```
