# Rank-based estimation of the linear regression coefficients

For a transformed Gaussian linear model, compute point estimates of the
regression coefficients. This approach uses the ranks of the data and
does not require the transformation, but must expand the sample size to
`n^2` and thus can be slow.

## Usage

``` r
rank_approx(y, X)
```

## Arguments

- y:

  `n x 1` response vector

- X:

  `n x p` matrix of predictors (should not include an intercept!)

## Value

the estimated linear coefficients

## Examples

``` r
# Simulate some data:
dat = simulate_tlm(n = 200, p = 10, g_type = 'step')

# Point estimates for the linear coefficients:
theta_hat = suppressWarnings(
  rank_approx(y = dat$y,
              X = dat$X[,-1]) # remove intercept
) # warnings occur from glm.fit (fitted probabilities 0 or 1)

# Check: correlation with true coefficients
cor(dat$beta_true[-1], # excluding the intercept
    theta_hat)
#> [1] 0.9132406
```
