
<!-- README.md is generated from README.Rmd. Please edit that file -->

# wildmeta

<!-- badges: start -->
<!-- badges: end -->

Typical methods to conduct meta-analysis—-pooling effect sizes or
analyzing moderating effects with meta-regression—-work under the
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
corrections for multiple-contrast hypothesis tests. Tipton & Pustejovsky
(2015) found that the HTZ test, which is an extension of the CR2
correction method with the Satterthwaite degrees of freedom, controlled
Type 1 error rate adequately even when the number of studies was small.
However, Joshi, Pustejovsky & Beretvas (2021) showed, through
simulations, that the HTZ test can be conservative. We examined another
method, cluster wild bootstrapping (CWB), that has been studied in the
econometrics literature but not in the meta-analytic context. The
results of the simulations from Joshi, Pustejovsky & Beretvas (2021)
showed that CWB adequately controlled for Type 1 error rate and had more
power than the HTZ test especially for multiple-contrast hypothesis
tests.

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
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               boot_F
#> 1 0.22863228, 0.13755777, 2.64087417, 0.18689552, 0.49963378, 0.71058282, 5.31290884, 5.54883243, 0.08164397, 0.71299038, 1.41893115, 1.37862321, 2.88319324, 0.20776536, 0.92085103, 2.46258905, 0.42688882, 1.65040593, 0.16793914, 0.36057311, 0.30580921, 2.10151464, 1.50768085, 0.63606900, 1.47704329, 0.01779652, 0.40321165, 0.61241492, 0.64319286, 0.56601766, 1.65598009, 0.40885110, 1.58376663, 1.84767951, 1.48587044, 0.50871317, 3.14110626, 1.25332355, 1.04294169, 0.07302695, 0.02535013, 0.13825696, 4.67609758, 0.05131970, 1.53069042, 1.72563102, 0.79204688, 0.04306586, 0.22158408, 4.84852200, 4.24340012, 1.50874825, 0.63851216, 1.83146738, 0.14005002, 0.23590387, 2.13008469, 0.77675930, 0.26747100, 0.68711400, 0.52857014, 3.30759904, 1.25475433, 0.15957752, 2.57873179, 0.83717477, 0.26854859, 0.04612052, 0.23126869, 0.20883892, 0.77200585, 1.29439561, 0.15466512, 1.26987033, 0.49652006, 5.64131939, 3.52993485, 2.02995561, 1.03987437, 2.50155770, 0.49929644, 0.22416879, 2.91943974, 0.15437753, 3.61217573, 0.53343871, 2.99561251, 0.32296307, 2.25466355, 0.14481443, 0.19671636, 0.45700715, 1.56756878, 0.23023245, 1.05863130, 0.11652668, 0.95194227, 0.52485875, 0.98264764
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

# References

Canty A. & Ripley B. (2020). boot: Bootstrap R (S-Plus) Functions. R
package version 1.3-25. <https://CRAN.R-project.org/package=boot>

DerSimonian, R., & Laird, N. (1983). Evaluating the effect of coaching
on SAT scores: A meta-analysis. Harvard Educational Review, 53(1), 1-1.

Fischer, A. & Roodman, D. (2021). fwildclusterboot: Fast Wild Cluster
Bootstrap Inference for Linear Regression Models. Available from
<https://cran.r-project.org/package=fwildclusterboot>.

Graham N., Arai M., & Hagströmer, B (2016). multiwayvcov: Multi-Way
Standard Error Clustering. R package version 1.2.3.
<https://CRAN.R-project.org/package=multiwayvcov>

Hedges, L. V., Tipton, E., & Johnson, M. C. (2010). Robust variance
estimation in meta-regression with dependent effect size estimates.
Research Synthesis Methods, 1(1), 39–65.
<https://doi.org/10.1002/jrsm.5>

Heyman, M. (2019). lmboot: Bootstrap in Linear Models. R package version
0.0.1. <https://CRAN.R-project.org/package=lmboot>

Joshi, M., Pustejovsky, J. E., & Beretvas, S. N. (2021). Cluster wild
bootstrapping to handle dependent effect sizes in meta-Analysis with
small number of studies. Working paper.

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
