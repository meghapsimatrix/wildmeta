
<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src="man/figures/wildmeta_hex.png" align="right" alt="" width="180" />

# wildmeta

<!-- badges: start -->

[![R-CMD-check](https://github.com/meghapsimatrix/wildmeta/workflows/R-CMD-check/badge.svg)](https://github.com/meghapsimatrix/wildmeta/actions)
[![Codecov test
coverage](https://codecov.io/gh/meghapsimatrix/wildmeta/branch/main/graph/badge.svg)](https://app.codecov.io/gh/meghapsimatrix/wildmeta?branch=main)
[![CRAN
status](https://www.r-pkg.org/badges/version/wildmeta)](https://CRAN.R-project.org/package=wildmeta)
[![](http://cranlogs.r-pkg.org/badges/grand-total/wildmeta)](https://CRAN.R-project.org/package=wildmeta)
[![](http://cranlogs.r-pkg.org/badges/last-month/wildmeta)](https://CRAN.R-project.org/package=wildmeta)
<!-- badges: end -->

Typical methods to conduct meta-analysis—pooling effect sizes or
analyzing moderating effects with meta-regression—work under the
assumption that the effect size estimates are independent. However,
primary studies often report multiple estimates of effect sizes.
Presence of multiple effect sizes leads to dependence as the estimates
within each study are likely correlated (e.g., because the same
participants provide multiple outcome scores). The increasingly popular
method to handle such dependence, robust variance estimation (RVE),
results in inflated Type 1 error rate when the number of studies is
small (Hedges, Tipton & Johnson, 2010; Tipton, 2015).

Tipton (2015) and Tipton & Pustejovsky (2015) examined several small
sample correction methods. Tipton (2015) recommended CR2 type correction
for RVE as well as the use of Satterthwaite degrees of freedom for
single coefficient tests. Tipton & Pustejovsky (2015) examined
corrections for multiple-contrast hypothesis tests. The authors found
that the HTZ test, which is an extension of the CR2 correction method
with the Satterthwaite degrees of freedom, controlled Type 1 error rate
adequately even when the number of studies was small. However, Joshi,
Pustejovsky & Beretvas (2022) showed, through simulations, that the HTZ
test can be conservative. The authors examined another method, cluster
wild bootstrapping (CWB), that has been studied in the econometrics
literature but not in the meta-analytic context. The results of the
simulations from Joshi, Pustejovsky & Beretvas (2021) showed that CWB
adequately controlled for Type 1 error rate and provided higher power
than the HTZ test, especially for multiple-contrast hypothesis tests.

The goal of this package is to provide applied meta-analytic researchers
a set of functions with which they can conduct single coefficient tests
or multiple-contrast hypothesis tests using cluster wild bootstrapping.

## Installation

Install the latest release from CRAN:

``` r
install.packages("wildmeta")
```

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("meghapsimatrix/wildmeta")
```

## Example

The following example uses the `SATCoaching` dataset from the
`clubSandwich` package (Pustejovksy, 2020), originally from DerSimonian
and Laird (1983). The standardized mean differences represent the
effects of SAT coaching on SAT verbal (SATV) and/or SAT math (SATM)
scores.

The code below runs cluster wild bootstrapping to test the
multiple-contrast hypothesis that the effect of coaching does not differ
based on study type. The `study_type` variable indicates whether groups
compared in primary studies were matched, randomized, or non-equivalent.
The meta-regression model also controls for hours of coaching provided
(`hrs`) and whether the students took math or verbal test (`test`).
Below, we run a zero-intercept meta-regression model.

Below, we use the `robumeta::robu()` function to fit the full model. The
functions in our package work with models fit using `robumeta::robu()`
(Fisher, Tipton, & Zhipeng, 2017) and `metafor::rma.mv()` (Viechtbauer,
2010).

``` r
library(wildmeta)
library(clubSandwich)
library(robumeta)

full_model <- robu(d ~ 0 + study_type + hrs + test,
                   studynum = study,
                   var.eff.size = V,
                   small = FALSE,
                   data = SATcoaching)

Wald_test_cwb(full_model = full_model,
              constraints = constrain_equal(1:3),
              R = 12,
              seed = 20201210)
#>   Test Adjustment CR_type Statistic  R     p_val
#> 1  CWB        CR0     CR0   Naive-F 12 0.4166667
```

# Related Work

We want to recognize other packages that provide algorithms to conduct
bootstrapping. The
[`fwildclusterboot`](https://s3alfisc.github.io/fwildclusterboot/index.html)
package runs the fast version of cluster wild bootstrapping proposed by
Roodman et al. (2019) and is based on the Stata’s boottest package
(Fischer & Roodman, 2021). The package runs bootstrapping for linear
regression models and fixed effects models. To the best of our
knowledge, meta-analytic models, weighted linear regression models, and
multiple contrast hypothesis tests are not supported by the package.

The [`multiwayvcov`](https://CRAN.R-project.org/package=multiwayvcov)
package implements cluster robust variance estimation as well as
different types of cluster bootstrapping, including pair, residual, and
cluster wild bootstrapping (Graham, Arai & Hagströmer, 2016). For
cluster wild bootstrapping, Rademacher weights, Mammen weights, and
standard normal weights are available. The functions in the package
return cluster robust variance-covariance matrices. Our package, on the
other hand, outputs p-values from bootstrap hypothesis tests. Further,
our package is particularly created for applied meta-analysis whereas
the `multiwayvcov` can be used with any regression analysis involving
clusters.

The [`boot`](https://CRAN.R-project.org/package=boot) package can be
used to generate bootstrap replicates (Canty & Ripley, 2020). The main
function in the package, `boot()` requires a dataset and a function that
calculates the estimate of interest. The function returns the bootstrap
estimate. The package also contains another function, `boot.ci()` that
can be used to calculate bootstrap confidence intervals. The package
allows parallel programming. However, the `boot` package does not
provide functionality for cluster wild bootstrapping.

The [`lmboot`](https://CRAN.R-project.org/package=lmboot) package
implements pair, residual and wild bootstrapping for linear models
(Heyman, 2019). The package does not provide Rademacher or Mammen
weights for wild bootstrapping. The output of the function is a
bootstrap sampling distribution. Further, the package does not work with
clustered data.

# Acknowledgments

Huge shout out to my extremely talented friend
[RAZ](https://ms-raz.com/) for creating our hex sticker within minutes.
We are also extremely thankful to [Wolfgang
Viechtbauer](https://wvbauer.com/doku.php/home) for helping us solve
[this
issue](https://stat.ethz.ch/pipermail/r-help/2021-November/472977.html)
we had with the `update()` function. We also thank Mikkel Vembye for
testing our package and giving us very helpful feedback.

# References

Canty A. & Ripley B. (2020). boot: Bootstrap R (S-Plus) Functions. R
package version 1.3-25. <https://CRAN.R-project.org/package=boot>

DerSimonian, R., & Laird, N. (1983). Evaluating the effect of coaching
on SAT scores: A meta-analysis. Harvard Educational Review, 53(1), 1-1.

Fischer, A. & Roodman, D. (2021). fwildclusterboot: Fast Wild Cluster
Bootstrap Inference for Linear Regression Models. Available from
<https://cran.r-project.org/package=fwildclusterboot>.

Fisher, Z., Tipton, E., & Zhipeng, H. (2017). robumeta: Robust variance
meta-regression. Retrieved from
<https://CRAN.R-project.org/package=robumeta>

Graham N., Arai M., & Hagströmer, B (2016). multiwayvcov: Multi-Way
Standard Error Clustering. R package version 1.2.3.
<https://CRAN.R-project.org/package=multiwayvcov>

Hedges, L. V., Tipton, E., & Johnson, M. C. (2010). Robust variance
estimation in meta-regression with dependent effect size estimates.
Research Synthesis Methods, 1(1), 39–65.

Heyman, M. (2019). lmboot: Bootstrap in Linear Models. R package version
0.0.1. <https://CRAN.R-project.org/package=lmboot>

Joshi, M., Pustejovsky, J. E., & Beretvas, S. N. (2022). Cluster wild
bootstrapping to handle dependent effect sizes in meta-analysis with
small number of studies. Research Synthesis Methods.
<https://doi.org/10.1002/jrsm.1554>

Pustejovsky, J. E. (2020). clubSandwich: Cluster-robust (sandwich)
variance estimators with small-sample corrections \[R package version
0.4.2\]. R package version 0.4.2.
<https://CRAN.R-project.org/package=clubSandwich>

Tipton, E. (2015). Small sample adjustments for robust variance
estimation with meta-regression. Psychological Methods, 20(3), 375–393.

Tipton, E., & Pustejovsky, J. E. (2015). Small-Sample Adjustments for
Tests of Moderators and Model Fit Using Robust Variance Estimation in
Meta-Regression. Journal of Educational and Behavioral Statistics (Vol.
40).

Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor
package. Journal of Statistical Software, 36(3), 1–48.
