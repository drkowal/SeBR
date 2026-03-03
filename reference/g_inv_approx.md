# Approximate inverse transformation

Compute the inverse function of a transformation `g` based on a grid
search.

## Usage

``` r
g_inv_approx(g, t_grid)
```

## Arguments

- g:

  the transformation function

- t_grid:

  grid of arguments at which to evaluate the transformation function

## Value

A function which can be used for evaluations of the (approximate)
inverse transformation function.
