# Inverse Box-Cox transformation

Evaluate the inverse Box-Cox transformation. Negative values are
permitted.

## Usage

``` r
g_inv_bc(s, lambda)
```

## Arguments

- s:

  argument(s) at which to evaluate the function

- lambda:

  Box-Cox parameter

## Value

The evaluation(s) of the inverse Box-Cox function at the given input(s)
`s`.

## Note

Special cases include the identity transformation (`lambda = 1`), the
square-root transformation (`lambda = 1/2`), and the log transformation
(`lambda = 0`).

## Examples

``` r
# (Inverse) log-transformation:
g_inv_bc(1:5, lambda = 0); exp(1:5)
#> [1]   2.718282   7.389056  20.085537  54.598150 148.413159
#> [1]   2.718282   7.389056  20.085537  54.598150 148.413159

# (Inverse) square-root transformation: note the shift and scaling
g_inv_bc(1:5, lambda = 1/2); (1:5)^2
#> [1]  2.25  4.00  6.25  9.00 12.25
#> [1]  1  4  9 16 25
```
