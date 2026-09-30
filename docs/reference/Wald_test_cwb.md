# Calculate p-values with cluster wild bootstrapping for meta-regression models.

Calculate p-values for single coefficient and multiple contrast
hypothesis tests using cluster wild bootstrapping.

## Usage

``` r
Wald_test_cwb(
  full_model,
  constraints,
  R,
  cluster = NULL,
  auxiliary_dist = "Rademacher",
  adjust = "CR0",
  type = "CR0",
  test = "Naive-F",
  seed = NULL,
  future_args = NULL
)
```

## Arguments

- full_model:

  Model fit using
  [`robumeta::robu()`](https://rdrr.io/pkg/robumeta/man/robu.html),
  [`metafor::rma.mv()`](https://wviechtb.github.io/metafor/reference/rma.mv.html),
  or
  [`metafor::rma.uni()`](https://wviechtb.github.io/metafor/reference/rma.uni.html)
  that includes the full set of moderators in the meta-regression model.

- constraints:

  A q X p constraint matrix be tested. Alternately, a function to create
  such a matrix, specified using
  [`clubSandwich::constrain_equal()`](http://jepusto.github.io/clubSandwich/reference/constraint_matrices.md)
  or
  [`clubSandwich::constrain_zero()`](http://jepusto.github.io/clubSandwich/reference/constraint_matrices.md).

- R:

  Number of bootstrap replications.

- cluster:

  Vector of identifiers indicating which observations belong to the same
  cluster. If `NULL` (the default), then the clustering variable will be
  inferred based on the structure of `full_mod`.

- auxiliary_dist:

  Character string indicating the auxiliary distribution to be used for
  cluster wild bootstrapping, with available options: "Rademacher",
  "Mammen", "Webb six", "uniform", "standard normal". The default is set
  to "Rademacher." We recommend the Rademacher distribution for models
  that have at least 10 clusters. For models with less than 10 clusters,
  we recommend the use of "Webb six" distribution.

- adjust:

  Character string specifying which small-sample adjustment should be
  used to multiply the residuals by. The available options are `"CRO"`,
  `"CR1"`, `"CR2"`, `"CR3"`, or `"CR4"`, with a default of `"CRO"`.

- type:

  Character string specifying which small-sample adjustment is used to
  calculate the Wald test statistic. The available options are `"CRO"`,
  `"CR1"`, `"CR2"`, `"CR3"`, or `"CR4"`, with a default of `"CRO"`.

- test:

  Character string specifying which (if any) small-sample adjustment is
  used in calculating the test statistic. Default is `"Naive-F"`, which
  does not make any small-sample adjustment.

- seed:

  Optional seed value to ensure reproducibility.

- future_args:

  Optional list of additional arguments passed to the `future_*()`
  functions used in calculating results across bootstrap replications.
  Ignored if the future.apply package is not available.

## Value

A `data.frame` containing the name of the test, the adjustment used for
the bootstrap process, the type of variance-covariance matrix used, the
type of test statistic, the number of bootstrap replicates, and the
bootstrapped p-value.

## Examples

``` r
library(clubSandwich)
library(robumeta)

model <- robu(d ~ 0 + study_type + hrs + test,
             studynum = study,
              var.eff.size = V,
              small = FALSE,
              data = SATcoaching)

C_mat <- constrain_equal(1:3, coefs = coef(model))

Wald_test_cwb(full_model = model,
              constraints = C_mat,
              R = 12)
#>   Test Adjustment CR_type Statistic  R p_val
#> 1  CWB        CR0     CR0   Naive-F 12   0.5

# Equivalent, using constrain_equal()
Wald_test_cwb(full_model = model,
              constraints = constrain_equal(1:3),
              R = 12)
#>   Test Adjustment CR_type Statistic  R     p_val
#> 1  CWB        CR0     CR0   Naive-F 12 0.5833333
```
