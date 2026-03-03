# Compute the latent data CDF

Assuming a Gaussian latent data distribution (given x), compute the CDF
on a grid of points

## Usage

``` r
Fz_fun(z, weights = NULL, mean_vec = NULL, sd_vec)
```

## Arguments

- z:

  vector of points at which the CDF of z is evaluated

- weights:

  `n`-dimensional vector of weights; if NULL, assume 1/n

- mean_vec:

  `n`-dimensional vector of means; if NULL, assume mean zero

- sd_vec:

  `n`-dimensional vector of standard deviations

## Value

CDF of z evaluated at `z`
