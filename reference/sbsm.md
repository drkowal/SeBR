# Semiparametric Bayesian spline model

Monte Carlo sampling for Bayesian spline regression with an unknown
(nonparametric) transformation. Cubic B-splines are used with a prior
that penalizes roughness.

## Usage

``` r
sbsm(
  y,
  x = NULL,
  x_test = NULL,
  psi = NULL,
  laplace_approx = TRUE,
  fixedX = (length(y) >= 500),
  approx_g = FALSE,
  nsave = 1000,
  ngrid = 100,
  verbose = TRUE
)
```

## Arguments

- y:

  `n x 1` response vector

- x:

  `n x 1` vector of observation points; if NULL, assume equally-spaced
  on \[0,1\]

- x_test:

  `n_test x 1` vector of testing points; default is `x`

- psi:

  prior variance (inverse smoothing parameter); if NULL, sample this
  parameter

- laplace_approx:

  logical; if TRUE, use a normal approximation to the posterior in the
  definition of the transformation; otherwise the prior is used

- fixedX:

  logical; if TRUE, treat the design as fixed (non-random) when sampling
  the transformation; otherwise treat covariates as random with an
  unknown distribution

- approx_g:

  logical; if TRUE, apply large-sample approximation for the
  transformation

- nsave:

  number of Monte Carlo simulations

- ngrid:

  number of grid points for inverse approximations

- verbose:

  logical; if TRUE, print time remaining

## Value

a list with the following elements:

- `coefficients` the posterior mean of the regression coefficients

- `fitted.values` the posterior predictive mean at the test points
  `x_test`

- `post_theta`: `nsave x p` samples from the posterior distribution of
  the regression coefficients

- `post_ypred`: `nsave x n_test` samples from the posterior predictive
  distribution at `x_test`

- `post_g`: `nsave` posterior samples of the transformation evaluated at
  the unique `y` values

- `model`: the model fit (here, `sbsm`)

as well as the arguments passed in.

## Details

This function provides fully Bayesian inference for a transformed spline
regression model using Monte Carlo (not MCMC) sampling. The
transformation is modeled as unknown and learned jointly with the
regression function (unless `approx_g = TRUE`, which then uses a point
approximation). This model applies for real-valued data, positive data,
and compactly-supported data (the support is automatically deduced from
the observed `y` values). The results are typically unchanged whether
`laplace_approx` is TRUE/FALSE; setting it to TRUE may reduce
sensitivity to the prior, while setting it to FALSE may speed up
computations for very large datasets. By default, `fixedX` is set to
FALSE for smaller datasets (`n < 500`) and TRUE for larger datasets
(`n >= 500`).

## Examples

``` r
# \donttest{
# Simulate some data:
n = 200 # sample size
x = sort(runif(n)) # observation points

# Transform a noisy, periodic function:
y = g_inv_bc(
  sin(2*pi*x) + sin(4*pi*x) + rnorm(n),
             lambda = .5) # Signed square-root transformation

# Fit the semiparametric Bayesian spline model:
fit = sbsm(y = y, x = x)
#> [1] "3 sec remaining"
#> [1] "2 sec remaining"
#> [1] "Total time: 4 seconds"
names(fit) # what is returned
#>  [1] "coefficients"   "fitted.values"  "post_theta"     "post_ypred"    
#>  [5] "post_g"         "post_psi"       "model"          "y"             
#>  [9] "X"              "sample_psi"     "laplace_approx" "fixedX"        
#> [13] "approx_g"      

# Note: this is Monte Carlo sampling...no need for MCMC diagnostics!

# Plot the model predictions (point and interval estimates):
pi_y = t(apply(fit$post_ypred, 2, quantile, c(0.05, .95))) # 90% PI
plot(x, y, type='n', ylim = range(pi_y,y),
     xlab = 'x', ylab = 'y', main = paste('Fitted values and prediction intervals'))
polygon(c(x, rev(x)),c(pi_y[,2], rev(pi_y[,1])),col='gray', border=NA)
lines(x, y, type='p') # observed points

lines(x, fitted(fit), lwd = 3) # fitted curve
#> Error in xy.coords(x, y): 'x' and 'y' lengths differ
# }
```
