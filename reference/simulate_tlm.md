# Simulate a transformed linear model

Generate training data (X, y) and testing data (X_test, y_test) for a
transformed linear model. The covariates are correlated Gaussian
variables. A user-specified proportion (`prop_sig`) of the regression
coefficients are nonozero (= 1) and the rest are zero. There are
multiple options for the transformation, which define the support of the
data (see below).

## Usage

``` r
simulate_tlm(
  n,
  p,
  g_type = "beta",
  n_test = 1000,
  heterosked = FALSE,
  lambda = 1,
  prop_sig = 0.5
)
```

## Arguments

- n:

  number of observations in the training data

- p:

  number of covariates

- g_type:

  type of transformation; must be one of `beta`, `step`, or `box-cox`

- n_test:

  number of observations in the testing data

- heterosked:

  logical; if TRUE, simulate the latent data with heteroskedasticity

- lambda:

  Box-Cox parameter (only applies for `g_type = 'box-cox'`)

- prop_sig:

  proportion of signals (nonzero coefficients)

## Value

a list with the following elements:

- `y`: the response variable in the training data

- `X`: the covariates in the training data

- `y_test`: the response variable in the testing data

- `X_test`: the covariates in the testing data

- `beta_true`: the true regression coefficients

- `g_true`: the true transformation, evaluated at y

## Details

The transformations vary in complexity and support for the observed
data, and include the following options: `beta` yields marginally
Beta(0.1, 0.5) data supported on \[0,1\]; `step` generates a
locally-linear inverse transformation and produces positive data; and
`box-cox` refers to the signed Box-Cox family indexed by `lambda`, which
generates real-valued data with examples including identity,
square-root, and log transformations.

## Note

The design matrices `X` and `X_test` do not include an intercept and
there is no intercept parameter in `beta_true`. The location/scale of
the data are not identified in general transformed regression models, so
recovering them is not a goal.

## Examples

``` r
# Simulate data:
dat = simulate_tlm(n = 100, p = 5, g_type = 'beta')
names(dat) # what is returned
#> [1] "y"         "X"         "y_test"    "X_test"    "beta_true" "g_true"   
hist(dat$y, breaks = 25) # marginal distribution

```
