## Resubmission

This is a resubmission. This version implements functionals for calculating bootstrap F statistics, fixes bugs in parallel processing setup that led to NA results for models estimated with rma.mv() and updates updated internals of plot.Wald_test_wildmeta to avoid use of ggplot2::aes_string(), which is now deprecated.
 


## Test environments

* local OS Big Sur, R 4.2.1
* ubuntu 20.04.3 LTS (on Github), R devel, release, oldrelease
* macOS-latest (on Github), R release
* windows-latest (on Github), R release
* win-builder (devel, release, oldrelease)
* r-hub:
  * Windows Server 2008 R2 SP1, R-devel, 32/64 bit
  * Ubuntu Linux 16.04 LTS, R-release, GCC
  * Fedora Linux, R-devel, clang, gfortran
  * Debian Linux, R-devel, GCC

## R CMD check results

There were no ERRORs or WARNINGs. There was one NOTE

* Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.1002/jrsm.1554
    From: inst/doc/cwbmeta.html
          README.md
    Status: 503
    Message: Service Unavailable
  URL: https://doi.org/10.1002/jrsm.5
    From: inst/doc/cwbmeta.html
    Status: 503
    Message: Service Unavailable
    
* Found the following (possibly) invalid DOIs:
  DOI: 10.1002/jrsm.1554
    From: DESCRIPTION
    Status: Service Unavailable
    Message: 503


The flagged URLs and DOI are correct.


## revdepcheck results

We checked 1 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

