# Compute the sum-squared-residuals term under Zellner's g-prior

These sum-squared-residuals (SSR) arise in the variance (or precision)
term under 1) Zellner's g-prior on the coefficients and a Gamma prior on
the error precision and 2) marginalization over the coefficients.

## Usage

``` r
SSR_gprior(y, X = NULL, psi)
```

## Arguments

- y:

  vector of response variables

- X:

  matrix of covariates; if NULL, return `sum(y^2)`

- psi:

  prior variance (g-prior)

## Value

a positive scalar
