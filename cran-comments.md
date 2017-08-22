## Test environments
* local ubuntu 16.04 install, R 3.4.1
* ubuntu 12.04 (on travis-ci), R 3.4.1
* local windows 10, R 3.4.1

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking R code for possible problems ... NOTE
  predict.ASDA: no visible global function definition for ‘predict’
  Undefined global functions or variables:
    predict
  Consider adding
    importFrom("stats", "predict")
  to your NAMESPACE file.
