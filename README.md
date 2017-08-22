
![alt tag](https://travis-ci.org/gumeo/accSDA.svg?branch=master)

# accSDA
## Accelerated Sparse Discriminant Analysis

This is the `R`-package accompanying the paper [Proximal Methods for Sparse Optimal Scoring and Discriminant Analysis](https://arxiv.org/pdf/1705.07194.pdf).

This package is currently under development, although most of the functionality is there already! You can now do sparse discriminant analysis with the package, but the visualization tools are being implemented and tested.

# Installation

You can install the package from CRAN or for the development version, you can install directly from github.

To install packages from github you need the `devtools` package. So install that if you haven't gotten it already!

Now you can proceed to install the package:
```R
library(devtools)
install_github("gumeo/accSDA")
library(accSDA)
```
And now you can start playing around with the package!

# Iris tutorial

The following is an example on how one could use the package on Fisher's Iris dataset. I choose the Iris dataset because most people are familiar with it. Other examples with *p>>n* examples will arive later!

```R
# Prepare training and test set
train <- c(1:40,51:90,101:140)
Xtrain <- iris[train,1:4]

# normalize is a function in the package
nX <- normalize(Xtrain)
Xtrain <- nX$Xc
Ytrain <- iris[train,5]
Xtest <- iris[-train,1:4]
Xtest <- normalizetest(Xtest,nX)
Ytest <- iris[-train,5]
     
# Define parameters for SDAD, i.e. ADMM optimization method
# Also try the SDAP and SDAAP methods, look at the documentation
# to read more about the parameters!
Om <- diag(4)+0.1*matrix(1,4,4) #elNet coef mat
gam <- 0.01
lam <- 0.01
method <- "SDAD"
q <- 2
control <- list(PGsteps = 100,
                PGtol = c(1e-5,1e-5),
                mu = 1,
                maxits = 100,
                tol = 1e-3,
                quiet = FALSE)
     
# Run the algorithm
res <- ASDA(Xt = Xtrain,
            Yt = Ytrain,
            Om = Om,
            gam = gam ,
            lam = lam,
            q = q,
            method = method,
            control = control)
     
# Can also just use the defaults:
# Default optimization method is SDAAP, accelerated proximal gradient.
resDef <- ASDA(Xtrain,Ytrain)
```
Now that you have gotten some results, you want to test the performance on the test set! What comes out of the `ASDA` function is an S3 object of class `ASDA` and there is a predict method in the package to predict the outcome of the classifier on new data!

```R
preds <- predict(res, newdata = Xtest)
```

# Future plans

Coming releases will include more plotting and printing functionality for the `ASDA` objects. A C++ backend is also in the pipeline along with some further extensions to handle different types of data.
