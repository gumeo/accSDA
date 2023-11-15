## Test environments
* R version 4.1.2 (2021-11-01), Platform: aarch64-apple-darwin20 (64-bit), Running under: macOS 14.0
* rhub

## R CMD check results

There were no ERRORs or WARNINGs.

## Notes
I got an e-mail from Brian Ripley about the ggthemes dependency. This update
removes that dependency. There was also a new warning:

* checking S3 generic/method consistency ... WARNING
  function(X, ...)
SZVD_kFold_cv.default:
  function(X, Y, folds, gams, beta, D, q, maxits, tol, ztol, feat,
           penalty, quiet)
SZVD:
  function(train, ...)

SZVD.default:
  function(train, gamma, D, penalty, scaling, tol, maxits, beta, quiet)

ZVD:
  function(A, ...)
ZVD.default:
  function(A, scaling, get_DVs)

SZVDcv:
  function(Atrain, ...)
SZVDcv.default:
  function(Atrain, Aval, k, num_gammas, g_mults, D, sparsity_pen,
           scaling, penalty, beta, tol, ztol, maxits, quiet)
SZVD_kFold_cv:
See section 'Generic functions and methods' in the 'Writing R
Extensions' manual.

I fixed this by adding the missing ellipsis to the default functions.
