# Compute all subsets of a set

Given a set of variables, compute the inclusion indicators for all
possible subsets.

## Usage

``` r
all_subsets(set)
```

## Arguments

- set:

  the set from which to compute all subsets (e.g., `1:p`)

## Value

a data frame where the rows indicate the `2^p` different subsets and the
columns indicate inclusion (logical) for each element in that subset

## References

Code adapted from
\<https://www.r-bloggers.com/2012/04/generating-all-subsets-of-a-set/\>
