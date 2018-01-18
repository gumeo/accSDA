#' @title Accelerated Sparse Discriminant Analysis
#'
#' @description Applies accelerated proximal gradient algorithm, proximal gradient algorithm
#' or alternating direction methods of multipliers algorithm to
#' the optimal scoring formulation of sparse discriminant analysis proposed
#' by Clemmensen et al. 2011.
#' \deqn{argmin{|(Y_t\theta-X_t\beta)|_2^2 + t|\beta|_1 + \lambda|\beta|_2^2}}{
#' argmin{|(Y_t*theta-X_t*b)|_2^2 + t*|beta|_1 + lambda*|beta|_2^2}}
#'
#' @param Xt n by p data matrix, (can also be a data.frame that can be coerced to a matrix)
#' @param Yt n by K matrix of indicator variables (Yij = 1 if i in class j).
#'     This will later be changed to handle factor variables as well.
#'     Each observation belongs in a single class, so for a given row/observation,
#'     only one element is 1 and the rest is 0.
#' @param Om p by p parameter matrix Omega in generalized elastic net penalty.
#' @param gam Regularization parameter for elastic net penalty.
#' @param lam Regularization parameter for l1 penalty, must be greater than zero.
#'     If cross-validation is used (\code{CV = TRUE}) then this must be a vector
#'     of length greater than one.
#' @param q Desired number of discriminant vectors.
#' @param control List of control arguments. See Details.
#' @param method This parameter selects which optimization method to use.
#'     It is specified as a character vector which can be one of the three values
#'     \describe{
#'       \item{\code{SDAP}}{Proximal gradient algorithm.}
#'       \item{\code{SDAAP}}{Accelerated proximal gradient algorithm.}
#'       \item{\code{SDAD}}{Alternating directions method of multipliers algorithm.}
#'     }
#'     Note that further parameters are passed to the function in the argument \code{control},
#'     which is a list with named components.
#' @param ... Additional arguments for \code{\link[MASS]{lda}} function in package MASS.
#'
#' @details The control list contains the following entries to further tune the
#'          algorithms.
#'          \describe{
#'            \item{\code{PGsteps}}{Maximum number if inner proximal gradient/ADMM
#'            algorithm for finding beta. Default value is 1000.}
#'            \item{\code{PGtol}}{Stopping tolerance for inner method. If the method is \code{SDAD},
#'            then this must be a vector of two values, absolute (first element) and relative
#'            tolerance (second element). Default value is 1e-5 for both absolute and
#'            relative tolerances.}
#'            \item{\code{maxits}}{Number of iterations to run. Default value is 250.}
#'            \item{\code{tol}}{Stopping tolerance. Default value is 1e-3.}
#'            \item{\code{mu}}{Penalty parameter for augmented Lagrangian term,
#'            must be greater than zero and only needs to be specified when using
#'            method \code{SDAD}. Default value is 1.}
#'            \item{\code{CV}}{Logical value which is \code{TRUE} if cross validation is supposed to be
#'            performed. If cross-validation is performed, then lam should be specified as
#'            a vector containing the regularization values to be tested. Default value is \code{FALSE}.}
#'            \item{\code{folds}}{Integer determining the number of folds in cross-validation. Not needed
#'            if CV is not specified. Default value is 5.}
#'            \item{\code{feat}}{Maximum fraction of nonzero features desired in validation scheme. Not needed
#'            if CV is not specified. Default value is 0.15.}
#'            \item{\code{quiet}}{Set to \code{FALSE} if status updates are supposed to be printed to the R console.
#'            Default value is \code{TRUE}. Note that this triggers a lot of printing to the console.}
#'            \item{\code{ordinal}}{Set to \code{TRUE} if the labels are ordinal. Only available for methods
#'            \code{SDAAP} and \code{SDAD}.}
#'            \item{\code{initTheta}}{Option to set the initial theta vector, by default it is a vector of all ones
#'            for the first theta.}
#'            \item{\code{bt}}{Logical indicating whether backtracking should be used, only applies to
#'            the Proximal Gradient based methods. By default, backtracking is not used.}
#'            \item{\code{L}}{Initial estimate for Lipshitz constant used for backtracking. Default value is 0.25.}
#'            \item{\code{eta}}{Scalar for Lipshitz constant. Default value is 1.25.}
#'          }
#'
#' @return \code{ASDA} returns an object of \code{\link{class}} "\code{ASDA}" including a list
#' with the following named components:
#'
#' \describe{
#'   \item{\code{call}}{The matched call.}
#'   \item{\code{B}}{p by q matrix of discriminant vectors, i.e. sparse loadings.}
#'   \item{\code{Q}}{K by q matrix of scoring vectors, i.e. optimal scores.}
#'   \item{\code{varNames}}{Names of the predictors used, i.e. column names of Xt.}
#'   \item{\code{origP}}{Number of variables in Xt.}
#'   \item{\code{fit}}{Output from function \code{\link[MASS]{lda}} on projected data.
#'   This is \code{NULL} the trivial solution is found, i.e. B is all zeroes. Use
#'   lower values of \code{lam} if that is the case.}
#'   \item{\code{classes}}{The classes in Yt.}
#'   \item{\code{lambda}}{The lambda/\code{lam} used, best value found by cross-
#'   validation if \code{CV} is \code{TRUE}.}
#' }
#' @seealso \code{\link{SDAAP}}, \code{\link{SDAP}} and \code{\link{SDAD}}
#' @note The input matrix Xt should be normalized, i.e. each column corresponding to
#'     a variable should have its mean subtracted and scaled to unit length. The functions
#'     \code{\link{normalize}} and \code{\link{normalizetest}} are supplied for this purpose in the package.
#' @examples
#'     set.seed(123)
#'     # Prepare training and test set
#'     train <- c(1:40,51:90,101:140)
#'     Xtrain <- iris[train,1:4]
#'     nX <- normalize(Xtrain)
#'     Xtrain <- nX$Xc
#'     Ytrain <- iris[train,5]
#'     Xtest <- iris[-train,1:4]
#'     Xtest <- normalizetest(Xtest,nX)
#'     Ytest <- iris[-train,5]
#'
#'     # Define parameters for Alternating Direction Method of Multipliers (SDAD)
#'     Om <- diag(4)+0.1*matrix(1,4,4) #elNet coef mat
#'     gam <- 0.0001
#'     lam <- 0.0001
#'     method <- "SDAD"
#'     q <- 2
#'     control <- list(PGsteps = 100,
#'                     PGtol = c(1e-5,1e-5),
#'                     mu = 1,
#'                     maxits = 100,
#'                     tol = 1e-3,
#'                     quiet = FALSE)
#'
#'     # Run the algorithm
#'     res <- ASDA(Xt = Xtrain,
#'                 Yt = Ytrain,
#'                 Om = Om,
#'                 gam = gam ,
#'                 lam = lam,
#'                 q = q,
#'                 method = method,
#'                 control = control)
#'
#'     # Can also just use the defaults, which is Accelerated Proximal Gradient (SDAAP):
#'     resDef <- ASDA(Xtrain,Ytrain)
#'
#'     # Some example on simulated data
#'     # Generate Gaussian data on three classes with plenty of redundant variables
#'
#'     # This example shows the basic steps on how to apply this to data, i.e.:
#'     #  1) Setup training data
#'     #  2) Normalize
#'     #  3) Train
#'     #  4) Predict
#'     #  5) Plot projected data
#'     #  6) Accuracy on test set
#'
#'     P <- 300 # Number of variables
#'     N <- 50 # Number of samples per class
#'
#'     # Mean for classes, they are zero everywhere except the first 3 coordinates
#'     m1 <- rep(0,P)
#'     m1[1] <- 3
#'
#'     m2 <- rep(0,P)
#'     m2[2] <- 3
#'
#'     m3 <- rep(0,P)
#'     m3[3] <- 3
#'
#'     # Sample dummy data
#'     Xtrain <- rbind(MASS::mvrnorm(n=N,mu = m1, Sigma = diag(P)),
#'                    MASS::mvrnorm(n=N,mu = m2, Sigma = diag(P)),
#'                    MASS::mvrnorm(n=N,mu = m3, Sigma = diag(P)))
#'
#'     Xtest <- rbind(MASS::mvrnorm(n=N,mu = m1, Sigma = diag(P)),
#'                    MASS::mvrnorm(n=N,mu = m2, Sigma = diag(P)),
#'                    MASS::mvrnorm(n=N,mu = m3, Sigma = diag(P)))
#'
#'     # Generate the labels
#'     Ytrain <- factor(rep(1:3,each=N))
#'     Ytest <- Ytrain
#'
#'     # Normalize the data
#'     Xt <- accSDA::normalize(Xtrain)
#'     Xtrain <- Xt$Xc # Use the centered and scaled data
#'     Xtest <- accSDA::normalizetest(Xtest,Xt)
#'
#'     # Train the classifier and increase the sparsity parameter from the default
#'     # so we penalize more for non-sparse solutions.
#'     res <- accSDA::ASDA(Xtrain,Ytrain,lam=0.01)
#'
#'     # Plot the projected training data, it is projected to
#'     # 2-dimension because we have 3 classes. The number of discriminant
#'     # vectors is maximum number of classes minus 1.
#'     XtrainProjected <- Xtrain%*%res$beta
#'
#'     plot(XtrainProjected[,1],XtrainProjected[,2],col=Ytrain)
#'
#'     # Predict on the test data
#'     preds <- predict(res, newdata = Xtest)
#'
#'     # Plot projected test data with predicted and correct labels
#'     XtestProjected <- Xtest%*%res$beta
#'
#'     plot(XtestProjected[,1],XtestProjected[,2],col=Ytest,
#'          main="Projected test data with original labels")
#'     plot(XtestProjected[,1],XtestProjected[,2],col=preds$class,
#'          main="Projected test data with predicted labels")
#'
#'     # Calculate accuracy
#'     sum(preds$class == Ytest)/(3*N) # We have N samples per class, so total 3*N
#'
#' @rdname ASDA
#' @export ASDA
ASDA <- function (Xt, ...) UseMethod("ASDA")


#' @return \code{NULL}
#'
#' @rdname ASDA
#' @method ASDA default
ASDA.default <- function(Xt, Yt, Om = diag(p), gam = 1e-3, lam = 1e-6, q = K-1, method = "SDAAP", control=list(), ...){
  # This is straight from nnet:::formula
  # This is used later to handle Yt as a factor
  class.ind <- function(cl) {
    Ik=diag(length(levels(cl)))
    x=Ik[as.numeric(cl),]
    dimnames(x) <- list(names(cl), levels(cl))
    x
  }
  if(missing(Yt)){
    stop("You must specify labels/classes Yt to use this function!")
  }
  if(is.factor(Yt))
  {
    classes <- levels(Yt)
    factorY <- Yt
    Yt <- class.ind(Yt)
  } else {
    classes <- colnames(Yt)
    if(is.null(classes)){
      # No column names, generate numbers instead...
      classes <- paste(1:ncol(Yt))
    }
    factorY <- factor(colnames(Yt)[apply(Yt, 1, which.max)])
  }
  # Define K for default inputs
  K <- length(classes)
  p <- dim(Xt)[2]

  # Similar to the use of the optim function in stats from base
  ## Defaults for the control variables:
  con <- list(PGsteps = 1000,
              PGtol = 1e-5,
              maxits = 250,
              tol = 1e-3,
              mu = NA,
              CV = FALSE,
              folds = 5,
              feat = 0.15,
              quiet = TRUE,
              ordinal = FALSE,
              initTheta = matrix(1:K,nrow=K,ncol=1),
              bt = FALSE,
              L = 0.25,
              eta = 1.5)
  nmsC <- names(con)
  if(method == "SDAD"){
    # Set special defaults for SDAD
    con$mu <- 1
    con$PGtol <- c(1e-5,1e-5)
  }
  # Overwrite with user supplied input!
  con[(namc <- names(control))] <- control
  # Warn user of input that cannot be used!
  if(length(noNms <- namc[!namc %in% nmsC]))
    warning("unknown names in control: ", paste(noNms,collapse=", "))

  if(nrow(con$initTheta) != K | is.matrix(con$initTheta) == FALSE){
    stop('The initTheta parameter supplied must be a matrix with one
         column and the number of rows equal to the number fo classes.')
  }

  if(con$eta < 0){
    stop('Backtracking multiplier must be positive!')
  }

  # Make sure that the method is acelerated proximal gradient if we have ordinal data
  # Also inform that the CV version has not yet been implemented
  if(con$ordinal == TRUE){
    if(method == "SDAP"){
      stop('Ordinal SDA is only implemented for APG and ADMM')
    }
    if(con$CV == TRUE){
      stop('A cross-validation functionality has not been implemented for Ordinal SDA!')
    }
  }

  # Go through the inputs to verify that they are correct.
  if(missing(Xt)){
    stop("You must specify a data matrix Xt to use this function!")
  }
  predNames <- colnames(Xt)

  if(missing(Om)){
    print("Om not specified, using default diag(p).")
  }
  if(dim(Om)[1] != dim(Xt)[2]){
    stop("Om must be a p by p matrix, where p is the number of variables/columns in Xt.")
  }
  if(missing(gam)){
    print("gam not specified, using default value 1e-3")
  }
  if(gam < 0){
    stop("gam must be positive!")
  }
  if(missing(lam)){
    print("lam not specified, using default value 1e-6")
  }
  for(i in 1:length(lam)){
    if(lam[i]<0){
      stop("All elements of lam must be positive!")
    }
  }
  if(missing(q)){
    print("q not specified, default is the max number K-1, where K is number of classes in Yt")
  }
  if(q > as.numeric(K)-1){
    stop("q is too high, at most K-1 variates allowed, where K is the number of classes!")
  }
  PGsteps <- con$PGsteps
  if(PGsteps < 0){
    stop("PGsteps must be a positive integer!")
  }
  PGtol <- con$PGtol
  if(length(PGtol) == 1){
    if(PGtol < 0){
      stop("PGtol must be positive!")
    }
  } else{
    if(PGtol[1] < 0 | PGtol[2] < 0){
      stop("PGtol values must be positive!")
    }
    if(method == "SDAP" | method == "SDAAP"){
      stop("PGtol must be a single numeric value when method is not SDAD!")
    }
  }

  maxits <- con$maxits
  if(maxits < 0){
    stop("maxits must be positive")
  }
  tol <- con$tol
  if(tol < 0){
    stop("tol must be positive!")
  }
  if(missing(method)){
    print("method not specified, using default SDAAP. You can choose from: SDAP, SDAAP or SDAD")
  }
  if(!(method == "SDAP" | method == "SDAAP" | method == "SDAD")){
    stop(paste(method, " is not a valid method! Use SDAP, SDAAP or SDAD"))
  }
  mu <- con$mu
  if(method == "SDAD" & mu < 0){
    stop("mu must be positive")
  }
  if(method == "SDAD" & length(PGtol)!=2){
    stop("When using SDAD you can specify PGtol as a vector with
         two components, absolute and relative tolerance. E.g. PGtol <- c(1e-5,1e-5)")
  }
  CV <- con$CV
  folds <- con$folds
  feat <- con$feat
  quiet <- con$quiet
  if(CV){
    if(length(lam)<2){
      stop("If you want to do cross validation, then specify lam as a vector with
           the parameters to test!")
    }
    if(missing(folds)){
      stop("Specify the number of folds for cross-validation!")
    }
    if(missing(feat)){
      stop("Specify the desired proportion of features!")
    }
    if(feat < 0 | feat > 1){
      stop("feat must be between 0 and 1!")
    }
  }
  if(!quiet){
    print("Preliminary check of input finished!")
    print(paste("Using method: ",method))
    if(CV){
      print("Running cross-validation! This may take a while!")
    }
  }
  ###
  # Input has been checked
  ###
  # Use the selected method
  ###
  if(!CV){
    if(method == "SDAP"){
      if(con$bt == FALSE){
        res <- SDAP(Xt=Xt, Yt=Yt, Om=Om, gam=gam, lam=lam, q=q,
                    PGsteps=PGsteps, PGtol=PGtol, maxits=maxits,
                    tol=tol, initTheta=con$initTheta)
      }else{
        res <- SDAP(Xt=Xt, Yt=Yt, Om=Om, gam=gam, lam=lam, q=q,
                    PGsteps=PGsteps, PGtol=PGtol, maxits=maxits,
                    tol=tol, initTheta=con$initTheta, bt=con$bt,
                    L=con$L, eta=con$eta)
      }
    } else if(method == "SDAAP"){
      selector <- rep(1,dim(Xt)[2])
      if(con$ordinal == TRUE){
        selector[(dim(Xt)[2]-K+2):dim(Xt)[2]] <- rep(0,length((dim(Xt)[2]-K+2):dim(Xt)[2]))
      }
      if(con$bt==FALSE){
        res <- SDAAP(Xt=Xt, Yt=Yt, Om=Om, gam=gam, lam=lam, q=q, PGsteps=PGsteps,
                     PGtol=PGtol, maxits=maxits, tol=tol, selector=selector,
                     initTheta=con$initTheta)
      }else{
        res <- SDAAP(Xt=Xt, Yt=Yt, Om=Om, gam=gam, lam=lam, q=q, PGsteps=PGsteps,
                     PGtol=PGtol, maxits=maxits, tol=tol, selector=selector,
                     initTheta=con$initTheta, bt=con$bt,
                     L=con$L, eta=con$eta)
      }
    } else{ # method is SDAD, input has been checked
      selector <- rep(1,dim(Xt)[2])
      if(con$ordinal == TRUE){
        selector[(dim(Xt)[2]-K+2):dim(Xt)[2]] <- rep(0,length((dim(Xt)[2]-K+2):dim(Xt)[2]))
      }
      res <- SDAD(Xt=Xt, Yt=Yt, Om=Om, gam=gam, lam=lam, mu=mu, q=q, PGsteps=PGsteps,
                  PGtol=PGtol, maxits=maxits, tol=tol, selector=selector,
                  initTheta=con$initTheta)
    }
  } else{
    if(method == "SDAP"){
      if(con$bt == FALSE){
        res <- SDAPcv(X=Xt, Y=Yt, folds=folds, Om=Om, gam=gam,
                      lam=lam, q=q, PGsteps=PGsteps, PGtol=PGtol,
                      maxits=maxits, tol=tol, feat=feat,
                      quiet=quiet, initTheta=con$initTheta)
      }else{
        res <- SDAPcv(X=Xt, Y=Yt, folds=folds, Om=Om, gam=gam,
                      lam=lam, q=q, PGsteps=PGsteps, PGtol=PGtol,
                      maxits=maxits, tol=tol, feat=feat,
                      quiet=quiet, initTheta=con$initTheta, bt=con$bt,
                      L=con$L, eta=con$eta)
      }
    } else if(method == "SDAAP"){
      if(con$bt == FALSE){
        res <- SDAAPcv(X=Xt, Y=Yt, folds=folds, Om=Om, gam=gam,
                       lam=lam, q=q, PGsteps=PGsteps, PGtol=PGtol,
                       maxits=maxits, tol=tol, feat=feat,
                       quiet=quiet, initTheta=con$initTheta)
      }else{
        res <- SDAAPcv(X=Xt, Y=Yt, folds=folds, Om=Om, gam=gam,
                       lam=lam, q=q, PGsteps=PGsteps, PGtol=PGtol,
                       maxits=maxits, tol=tol, feat=feat,
                       quiet=quiet, initTheta=con$initTheta, bt=con$bt,
                       L=con$L, eta=con$eta)
      }
    } else{ # method is SDAD
      res <- SDADcv(X=Xt, Y=Yt, folds=folds, Om=Om, gam=gam,
                    lam=lam, mu=mu, q=q, PGsteps=PGsteps, PGtol=PGtol,
                    maxits=maxits, tol=tol, feat=feat,
                    quiet=quiet, initTheta=con$initTheta)
    }
  }
  ###
  # Done running the optimization
  ###
  # Use lda from MASS package on projected data
  ###
  if(!quiet){
    print("Algorithm has finished running, running lda on projected data and wrapping up!")
  }
  B <- res$B
  Q <- res$Q

  if(CV){
    lbest <- res$lbest
    lambest <- res$lambest
  } else{
    lambest <- lam
  }

  # Check if we get the trivial solution
  if (all(B==0)){
    # lda will throw error in this case
    # so we give a warning...
    warning('no non-zero elements - try other regularization parameters')
    lobj <- NULL
  }else{
    # Projected data
    sl <- Xt%*%B
    colnames(sl) <- paste("score", 1:ncol(sl), sep = "")
    lobj <- MASS::lda(sl, factorY, ...)
  }
  # Return structure
  structure(
    list(call = match.call(),
         beta = B,
         theta = Q,
         varNames = predNames,
         origP = dim(Xt)[2],
         fit = lobj,
         classes = classes,
         lambda = lambest,
         CV = CV),
    class = "ASDA")
}

#' @export
ASDA.data.frame <- function(Xt, ...){
  res <- ASDA.default(Xt = structure(data.matrix(Xt), class = "matrix"), ...)
  cl <- match.call()
  cl[[1L]] <- as.name("ASDA")
  res$call <- cl
  res
}

#' @export
ASDA.matrix <- function(Xt, ...){
  res <- ASDA.default(Xt, ...)
  cl <- match.call()
  cl[[1L]] <- as.name("ASDA")
  res$call <- cl
  res
}

#' Predict method for sparse discriminant analysis
#'
#' Predicted values based on fit from the function \code{\link{ASDA}}. This
#' function is used to classify new observations based on their explanatory variables/features.
#'
#' @param object Object of class ASDA. This object is returned from the function \code{\link{ASDA}}.
#' @param newdata A matrix of new observations to classify.
#' @param ... Arguments passed to \code{\link[MASS]{predict.lda}}.
#' @return A list with components:
#'
#' \describe{
#'   \item{\code{class}}{The classification (a factor)}
#'   \item{\code{posterior}}{posterior probabilities for the classes}
#'   \item{\code{x}}{the scores}
#' }
#' @seealso \code{\link{SDAAP}}, \code{\link{SDAP}} and \code{\link{SDAD}}
#' @note The input matrix newdata should be normalized w.r.t. the normalization
#'     of the training data
#' @examples
#'     # Prepare training and test set
#'     train <- c(1:40,51:90,101:140)
#'     Xtrain <- iris[train,1:4]
#'     nX <- normalize(Xtrain)
#'     Xtrain <- nX$Xc
#'     Ytrain <- iris[train,5]
#'     Xtest <- iris[-train,1:4]
#'     Xtest <- normalizetest(Xtest,nX)
#'     Ytest <- iris[-train,5]
#'
#'     # Define parameters for SDAD
#'     Om <- diag(4)+0.1*matrix(1,4,4) #elNet coef mat
#'     gam <- 0.01
#'     lam <- 0.01
#'     method <- "SDAD"
#'     q <- 2
#'     control <- list(PGsteps = 100,
#'                     PGtol = c(1e-5,1e-5),
#'                     mu = 1,
#'                     maxits = 100,
#'                     tol = 1e-3,
#'                     quiet = FALSE)
#'
#'     # Run the algorithm
#'     res <- ASDA(Xt = Xtrain,
#'                 Yt = Ytrain,
#'                 Om = Om,
#'                 gam = gam ,
#'                 lam = lam,
#'                 q = q,
#'                 method = method,
#'                 control = control)
#'
#'     # Do the predictions on the test set
#'     preds <- predict(object = res, newdata = Xtest)
#' @rdname predict.ASDA
#' @export
predict.ASDA <- function(object, newdata = NULL, ...)
{
  if(!is.matrix(newdata)) newdata <- as.matrix(newdata)
  if(!is.null(object$varNames) & length(object$varNames) > 0)
  {
    newdata <- newdata[, object$varNames, drop = FALSE]
  } else {
    if(ncol(newdata) != object$origP) stop("dimensions of training and testing X different")
    #newdata <- newdata[, object$varIndex, drop = FALSE]

  }
  xnew <- newdata %*%object$beta
  pred <- stats::predict(object$fit,xnew, ...)
  pred
}

#' Print method for ASDA object
#'
#' Prints a summary of the output from the \code{\link{ASDA}} function. The
#' output summarizes the discriminant analysis in human readable format.
#'
#' @param x Object of class ASDA. This object is returned from the function \code{\link{ASDA}}.
#' @param digits Number of digits to show in printed numbers.
#' @param numshow Number of best ranked variables w.r.t. to their absolute coefficients.
#' @param ... arguments passed to or from other methods.
#'
#' @return An invisible copy of \code{x}.
#' @seealso \code{\link{ASDA}}, \code{\link{predict.ASDA}} and \code{\link{SDAD}}
#' @examples
#'     # Prepare training and test set
#'     train <- c(1:40,51:90,101:140)
#'     Xtrain <- iris[train,1:4]
#'     nX <- normalize(Xtrain)
#'     Xtrain <- nX$Xc
#'     Ytrain <- iris[train,5]
#'     Xtest <- iris[-train,1:4]
#'     Xtest <- normalizetest(Xtest,nX)
#'     Ytest <- iris[-train,5]
#'
#'     # Run the algorithm
#'     resDef <- ASDA(Xtrain,Ytrain)
#'
#'     # Print
#'     print(resDef)
#' @rdname print.ASDA
#' @export
print.ASDA <- function(x, digits = max(3, getOption("digits") - 3), numshow = 5, ...)
{
  # Print the original call
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  # Gather classes into a string
  classInfo <- paste(x$classes, collapse = ", ")

  cat("lambda =", format(x$lambda, digits = digits),
      "\nclasses =", classInfo,
      "\n\n")

  top <- if(!is.null(x$varNames)) x$varNames else paste("Predictor", 1:dim(x$beta)[1], sep = "")
  varOrder <- if(is.matrix(x$beta)) order(apply(abs(x$beta), 1, sum)) else order(abs(x$beta))
  top <- top[varOrder]
  top <- top[1:min(numshow, length(top))]
  top <- paste(top, collapse = ", ")
  topStr <- paste("Top", min(numshow, dim(x$beta)[1]),"predictors:\t")

  cat(topStr,
      top,
      "\n",
      sep = "")

  invisible(x)
}

# TODO:
#-----------------------------------------------------------------------------------------------
#       Implement barplot
#-----------------------------------------------------------------------------------------------
#       Implement plot, also status diagnostic plot w.r.t. classification rates
#-----------------------------------------------------------------------------------------------
#       Implement summary
#       Likely just call print to start with
#-----------------------------------------------------------------------------------------------
#       Add verbose for intermediate printing of function value, implement this!
#-----------------------------------------------------------------------------------------------
# This is maybe not relevant, better for shorter tests
#       Test the package microbenchmark to measure times.
#-----------------------------------------------------------------------------------------------
#       Include Penicillin data, document it and get column names from Line.
#-----------------------------------------------------------------------------------------------
#       Add option to do CV in parallel
#-----------------------------------------------------------------------------------------------
# Find good profiler for this
#       Maybe do profiling to spot potential for C++ implementations of some parts.
#-----------------------------------------------------------------------------------------------
#       Do the same tests as Summer.
#-----------------------------------------------------------------------------------------------
#       Create an onload message, concerning maybe citation and other stuff
#-----------------------------------------------------------------------------------------------
# Finished!
#       Make input for ASDA more user friendly.
#          In this regard look at control in:
#          https://github.com/lgautier/R-3-0-branch-alt/blob/master/src/library/stats/R/optim.R
#          There are several things that can easily be added to
#          a control list variable.
#       Done: Needs testing! Done, works!
#-----------------------------------------------------------------------------------------------
# Finished! Maybe add the print from MASS::print.lda as well
#       Implement print
#       Add proportion of nonzero components in betas
#-----------------------------------------------------------------------------------------------
# Finished for now, need to revisit this for other functions.
#       Add examples into documentation
#-----------------------------------------------------------------------------------------------
# Remember to move this list to somewhere on github and create news.md and other files
#-----------------------------------------------------------------------------------------------
