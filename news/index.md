# Changelog

## SeBR 1.1.0

CRAN release: 2025-06-16

- [`bb()`](https://drkowal.github.io/SeBR/reference/bb.md) and
  [`hbb()`](https://drkowal.github.io/SeBR/reference/hbb.md) now omit
  scaling by n/(n+1). The scaling is part of the broader `SeBR` methods
  but is not part of the Bayesian bootstrap or hierarchical Bayesian
  bootstrap.
