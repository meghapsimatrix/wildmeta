
<!-- README.md is generated from README.Rmd. Please edit that file -->

# wildmeta

<!-- badges: start -->
<!-- badges: end -->

Primary studies in education and social sciences often contain multiple
effect sizes (Hedges, Tipton & Johnson, 2010). Presence of multiple
effect sizes leads to dependence as the estimates within each study are
likely correlated (e.g., because the same participants provide outcome
scores). The increasingly popular method to handle such dependence,
robust variance estimation (RVE), results in inflated Type 1 error rate
when the number of studies is small (Tipton, 2015). Tipton (2015) and
Tipton and Pustejovsky (2015) have recommended CR2 type correction to
RVE with Satterthwaite degrees of freedom. Tipton (2015) and Tipton and
Pustejovsky (2015) showed that the CR2 correction controls Type 1 error
rate adequately even when the number of studies is small.

Through simulations that I ran for my dissertation, I showed the the CR2
test with Satterthwaite degrees of freedom can be conservative. I
examined another method, cluster wild bootstrapping (CWB), that has been
studied in the econometrics literature but not in the meta-analytic
context. The results of my dissertation simulations showed that CWB
adequately controls for Type 1 error rate and has more power than the
CR2 test with Satterthwaite degrees of freedom.

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

This is a basic example which shows you how to solve a common problem:

``` r
library(wildmeta)
## basic example code
```

# References

Hedges, L. V., Tipton, E., & Johnson, M. C. (2010). Robust variance
estimation in meta-regression with dependent effect size estimates.
Research Synthesis Methods, 1(1), 39–65.
<https://doi.org/10.1002/jrsm.5>

Tipton, E. (2015). Small sample adjustments for robust variance
estimation with meta-regression. Psychological Methods, 20(3), 375–393.
<https://doi.org/10.1037/met0000011>

Tipton, E., & Pustejovsky, J. E. (2015). Small-Sample Adjustments for
Tests of Moderators and Model Fit Using Robust Variance Estimation in
Meta-Regression. Journal of Educational and Behavioral Statistics (Vol.
40). <https://doi.org/10.3102/1076998615606099>
