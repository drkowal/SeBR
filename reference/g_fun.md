# Compute the transformation

Given the CDFs of z and y, compute a smoothed function to evaluate the
transformation

## Usage

``` r
g_fun(y, Fy_eval, z, Fz_eval)
```

## Arguments

- y:

  vector of points at which the CDF of y is evaluated

- Fy_eval:

  CDF of y evaluated at `y`

- z:

  vector of points at which the CDF of z is evaluated

- Fz_eval:

  CDF of z evaluated at `z`

## Value

A smooth monotone function which can be used for evaluations of the
transformation.
