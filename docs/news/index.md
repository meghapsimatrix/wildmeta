# Changelog

## wildmeta 0.3.2

CRAN release: 2023-03-08

- Fixed argument names inconsistency

## wildmeta 0.3.1

CRAN release: 2023-02-15

- Implemented functionals for calculating bootstrap F statistics.
- Fixed bugs in parallel processing setup that led to NA results for
  models estimated with rma.mv().
- Updated internals of plot.Wald_test_wildmeta to avoid use of
  ggplot2::aes_string(), which is now deprecated.

## wildmeta 0.3.0

CRAN release: 2022-10-07

- Added methods for rma.uni objects.

## wildmeta 0.2.0

CRAN release: 2022-06-29

- Added parallel processing via the future package (using future.apply).

## wildmeta 0.1.0

CRAN release: 2022-02-02

- First version
