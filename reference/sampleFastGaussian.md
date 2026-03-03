# Sample a Gaussian vector using Bhattacharya et al. (2016)

Sample from N(mu, Sigma) where Sigma = solve(crossprod(Phi) + solve(D))
and mu = Sigma\*crossprod(Phi, alpha):

## Usage

``` r
sampleFastGaussian(Phi, Ddiag, alpha)
```

## Arguments

- Phi:

  `n x p` matrix (of predictors)

- Ddiag:

  `p x 1` vector of diagonal components (of prior variance)

- alpha:

  `n x 1` vector (of data, scaled by variance)

## Value

Draw from N(mu, Sigma), which is `p x 1`, and is computed in `O(n^2*p)`

## Note

Assumes D is diagonal, but extensions are available

## References

Bhattacharya, Chakraborty, and Mallick (2016,
\<https://doi.org/10.1093/biomet/asw042\>)
