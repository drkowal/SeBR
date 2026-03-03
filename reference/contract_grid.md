# Grid contraction

Contract the grid if the evaluation points exceed some threshold. This
removes the corresponding z values. We can add points back to achieve
the same (approximate) length.

## Usage

``` r
contract_grid(z, Fz, lower, upper, add_back = TRUE, monotone = TRUE)
```

## Arguments

- z:

  grid points (ordered)

- Fz:

  function evaluated at those grid points

- lower:

  lower threshold at which to check Fz

- upper:

  upper threshold at which to check Fz

- add_back:

  logical; if true, expand the grid to (about) the original size

- monotone:

  logical; if true, enforce monotonicity on the expanded grid

## Value

a list containing the grid points `z` and the (interpolated) function
`Fz` at those points
