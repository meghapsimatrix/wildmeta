## Resubmission
This is a resubmission. This is a minor release which adds parallel processing capabilities to run cluster wild bootstrapping.


## Test environments

* local OS Big Sur, R 4.1.2
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

There were no ERRORs or WARNINGs or NOTEs. 


## rhub and win_devel check results

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
