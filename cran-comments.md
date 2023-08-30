## Resubmission 30/08/2023 v.1.1.1

* Fix notes about long running examples by wrapping
  in \donttest{} blocks.
* Previous submission (v.1.1.0) also had the following
  notes from `r-devel-linux-x86_64-debian-gcc` flavor
  check:
    - Running R code in 'testthat.R' had CPU time
      7.6 times elapsed time
    - Re-building vignettes had CPU time
      10.8 times elapsed time
  I could not replicate these in any other environment:
    - locally using R CMD Check (Windows),
    - on r_hub (Windows, macos, and linux)
    - using R CMD Check in GitHub Actions runners
      (Windows, Ubuntu, macos)
    - using `devtools::check_win_devel()` and
      `check_macos_release()`
  I believe these to be false positives likely due to
  hardware or current load at the time of the check
  in case they are not, please provide advice how to fix.

## R CMD check results

0 errors | 0 warnings | 1 note

* Found possibly invalid URLs -- the URLs were
  checked to be valid and working, seems to be
  a false positive


## Reverse dependencies

* `survex` does not have any reverse dependencies
