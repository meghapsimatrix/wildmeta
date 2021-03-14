
<!-- README.md is generated from README.Rmd. Please edit that file -->

# wildmeta

<!-- badges: start -->
<!-- badges: end -->

Primary studies in education and social sciences often contain multiple
effect sizes (Hedges, Tipton & Johnson, 2010). Presence of multiple
effect sizes leads to dependence as the estimates within each study are
likely correlated (e.g., because the same participants provide multiple
outcome scores). The increasingly popular method to handle such
dependence, robust variance estimation (RVE), results in inflated Type 1
error rate when the number of studies is small (Hedges, Tipton &
Johnson, 2010; Tipton, 2015).

Tipton (2015) and Tipton and Pustejovsky (2015) recommended a small
sample correction, the HTZ test (CR2 correction with Satterthwaite
degrees of freedom). The HTZ test has been shown to control Type 1 error
rate adequately even when the number of studies is small (Tipton, 2015;
Tipton & Pustejovsky, 2015). Through simulations that I ran for my
dissertation, I showed the the HTZ test can be conservative. I examined
another method, cluster wild bootstrapping (CWB), that has been studied
in the econometrics literature but not in the meta-analytic context. The
results of my dissertation simulations showed that CWB adequately
controls for Type 1 error rate and has more power than the HTZ test.

The goal of this package is to provide applied meta-analytic researchers
a function with which they can conduct single coefficient tests or
multiple-contrast hypothesis tests using cluster wild bootstrapping.

## Installation

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
scores. The data contains the `study_type` variable indicating whether
groups compared in primary studies were matched, randomized, or
non-equivalent. Below, we use the `robu()` function from the `robumeta`
package to fit the full model. The we run cluster wild bootstrapping to
test the multiple-contrast hypothesis that the effect of coaching does
not differ based on study type.

``` r
library(wildmeta)
library(clubSandwich)
library(robumeta)

set.seed(12102020)


full <- robu(d ~ study_type,
             studynum = study,
             var.eff.size = V,
             small = FALSE,
             data = SATcoaching)


cwb(full_model = full,
    indices = 2:3,
    R = 99)
#>   test working_model     p_val
#> 1  CWB            CE 0.5252525
```

# Related Work

We want to recognize other packages that provide algorithms to conduct
bootstrapping. The `multiwayvcov` package implements cluster robust
variance estimation as well as different types of cluster bootstrapping,
including pair, residual, and cluster wild bootstrapping (Graham, Arai &
Hagströmer, 2016). For cluster wild bootstrapping, Rademacher weights,
Mammen weights, and standard normal weights are available. The functions
in the package return cluster robust variance-covariance matrices. Our
package, on the other hand, outputs p-values from bootstrap hypothesis
tests. Further, our package is particularly created for applied
meta-analysis whereas the `multiwayvcov` can be used with any regression
analysis involving clusters.

The `boot` package can be used to generate bootstrap replicates (Canty &
Ripley, 2020). The main function in the package, `boot()` requires a
dataset and a function that calculates the estimate of interest. The
function returns the bootstrap estimate. The package also contains
another function, `boot.ci()` that can be used to calculate bootstrap
confidence intervals. The package allows parallel programming. However,
the `boot` package does not provide functionality for cluster wild
bootstrapping.

The `lmboot` package implements pair, residual and wild bootstrapping
for linear models (Heyman, 2019). The package does not provide
Rademacher or Mammen weights for wild bootstrapping. The output of the
function is a bootstrap sampling distribution. Further, the package does
not work with clustered data.

# References

Canty A. & Ripley B. (2020). boot: Bootstrap R (S-Plus) Functions. R
package version 1.3-25. <https://CRAN.R-project.org/package=boot>

DerSimonian, R., & Laird, N. (1983). Evaluating the effect of coaching
on SAT scores: A meta-analysis. Harvard Educational Review, 53(1), 1-1

Hedges, L. V., Tipton, E., & Johnson, M. C. (2010). Robust variance
estimation in meta-regression with dependent effect size estimates.
Research Synthesis Methods, 1(1), 39–65.
<https://doi.org/10.1002/jrsm.5>

Graham N., Arai M., & Hagströmer, B (2016). multiwayvcov: Multi-Way
Standard Error Clustering. R package version 1.2.3.
<https://CRAN.R-project.org/package=multiwayvcov>

Heyman, M. (2019). lmboot: Bootstrap in Linear Models. R package version
0.0.1. <https://CRAN.R-project.org/package=lmboot>

Pustejovsky, J. E. (2020). clubSandwich: Cluster-robust (sandwich)
variance estimators with small-sample corrections \[R package version
0.4.2\]. R package version 0.4.2.
<https://CRAN.R-project.org/package=clubSandwich>

Tipton, E. (2015). Small sample adjustments for robust variance
estimation with meta-regression. Psychological Methods, 20(3), 375–393.
<https://doi.org/10.1037/met0000011>

Tipton, E., & Pustejovsky, J. E. (2015). Small-Sample Adjustments for
Tests of Moderators and Model Fit Using Robust Variance Estimation in
Meta-Regression. Journal of Educational and Behavioral Statistics (Vol.
40). <https://doi.org/10.3102/1076998615606099>
