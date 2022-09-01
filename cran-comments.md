## Test environments
* local ubuntu 20.04 install, R 4.1.3
* rhub 
* winbuilder R Under development (unstable) (2022-08-31 r82783 ucrt)

## R CMD check results

There were no ERRORs or WARNINGs. 

There is one NOTE that is only found on Windows (Server 2022, R-devel 64-bit): 

```
* checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'
```
As noted in [R-hub issue #503](https://github.com/r-hub/rhub/issues/503), this could be due to a bug/crash in MiKTeX and can likely be ignored.

## Notes
I got an e-mail from Brian Ripley about an error in the current version on CRAN (1.1.1): https://cran.r-project.org/web/checks/check_results_accSDA.html

The error came from a LaTeX error in some internal documentation, this has now
been fixed.
