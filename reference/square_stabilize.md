# Numerically stabilize the squared elements

Given a vector to be squared, add a numeric buffer for the elements very
close to zero.

## Usage

``` r
square_stabilize(vec)
```

## Arguments

- vec:

  vector of inputs to be squared

## Value

a vector of the same length as \`vec\`
