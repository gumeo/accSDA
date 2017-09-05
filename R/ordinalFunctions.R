#' @title Ordinal Accelerated Sparse Discriminant Analysis
#'
#' @description Applies accelerated proximal gradient algorithm to
#' the optimal scoring formulation of sparse discriminant analysis proposed
#' by Clemmensen et al. 2011. The problem is further casted to a binary
#' classification problem as described in "Learning to Classify Ordinal Data:
#' The Data Replication Method" by Cardoso and da Costa to handle the ordinal labels.
#' This function serves as a wrapper for the \code{\link{ASDA}} function, where the
#' appropriate data augmentation is performed. Since the problem is casted into
#' a binary classication problem, only a single discriminant vector comes from the
#' result. The first *p* entries correspond to the variables/coefficients for
#' the predictors, while the following K-1 entries correspond to biases for the
#' found hyperplane, to separate the classes. The resulting object is of class ordASDA
#' and has an accompanying predict function. The paper by Cardoso and dat Costa can
#' be found here: (http://www.jmlr.org/papers/volume8/cardoso07a/cardoso07a.pdf).
#'
#' @param Xt n by p data matrix, (can also be a data.frame that can be coerced to a matrix)
#' @param Yt vector of length n, equal to the number of samples. The classes should be
#'           1,2,...,K where K is the number of classes. Yt needs to be a numeric vector.
#' @param Om p by p parameter matrix Omega in generalized elastic net penalty, where
#'           p is the number of variables.
#' @param s We need to find a hyperplane that separates all classes with different biases.
#'          For each new bias we define a binary classification problem, where a maximum of
#'          s ordinal classes or contained in each of the two classes. A higher value of s means
#'          that more data will be copied in the data augmentation step. BY default s is 1.
#' @param gam Regularization parameter for elastic net penalty, must be greater than zero.
#' @param lam Regularization parameter for l1 penalty, must be greater than zero.
#' @param control List of control arguments further passed to ASDA. See \code{\link{ASDA}}.
#' @param ... Additional arguments for \code{\link{ASDA}} and \code{\link[MASS]{lda}}
#'            function in package MASS.
#'
#'
#' @return \code{ordASDA} returns an object of \code{\link{class}} "\code{ordASDA}" including a list
#' with the same components as an ASDA objects and:
#'
#' \describe{
#'   \item{\code{XN}}{Normalized data, used in the predict function to normalize test data.}
#'   \item{\code{h}}{Scalar value for biases.}
#'   \item{\code{K}}{Number of classes.}
#' }
#' @seealso \code{\link{ASDA}}.
#' @note The input matrix Xt is normalized in the function, so no prior normalization is needed.
#'      The functions \code{\link{normalize}} and \code{\link{normalizetest}} are used
#'      for this purpose and are supplied in the package.
#' @examples
#'     set.seed(123)
#'
#'     # You can play around with these values to generate some 2D data to test one
#'     numClasses <- 15
#'     sigma <- matrix(c(1,-0.2,-0.2,1),2,2)
#'     mu <- c(0,0)
#'     numObsPerClass <- 5
#'
#'     # Generate the data, can access with train$X and train$Y
#'     train <- accSDA::genDat(numClasses,numObsPerClass,mu,sigma)
#'     test <- accSDA::genDat(numClasses,numObsPerClass*2,mu,sigma)
#'
#'     # Visualize it, only using the first variable gives very good separation
#'     plot(train$X[,1],train$X[,2],col = factor(train$Y),asp=1,main="Training Data")
#'
#'     # Train the ordinal based model
#'     res <- accSDA::ordASDA(train$X,train$Y,s=2,h=1, gam=1e-6, lam=1e-3)
#'     vals <- predict(object = res,newdata = test$X) # Takes a while to run ~ 10 seconds
#'     sum(vals==test$Y)/length(vals) # Get accuracy on test set
#'     plot(test$X[,1],test$X[,2],col = factor(test$Y),asp=1,
#'          main="Test Data with correct labels")
#'     plot(test$X[,1],test$X[,2],col = factor(vals),asp=1,
#'          main="Test Data with predictions from ordinal classifier")
#'
#' @rdname ordASDA
#' @export ordASDA
ordASDA <- function (Xt, ...) UseMethod("ordASDA")


#' @return \code{NULL}
#'
#' @rdname ordASDA
#' @method ordASDA default
ordASDA.default <- function(Xt, Yt, s=1, Om, gam = 1e-3, lam = 1e-6, control,...){
  h <- 1
  if(missing(Yt)){
    stop('We need the ordinal labels Yt to run this function!')
  }
  K <- length(unique(Yt))
  if(K==2){
    stop("Only two types of labels, just use a binary classifier!")
  }
  if(missing(Om)){
    Om <- matrix(0,dim(Xt)[2]+K-1,dim(Xt)[2]+K-1)
    for(i in 1:(dim(Xt)[2])){
      Om[i,i] <- 1
    }
  }else{
    if(dim(Om)[1] != dim(Xt)[2]){
      stop('Om must be of dimension p by p, where p is the number of predictors!')
    }
    regMat <- matrix(0,dim(Xt)[2]+K-1,dim(Xt)[2]+K-1)
    regMat[1:dim(Xt)[2],1:dim(Xt)[2]] <- Om
    Om <- regMat
  }

  if(missing(control)){
    control <- list(ordinal=TRUE)
  }else{
    control$ordinal <- TRUE
  }

  # Normalize the data
  XN <- accSDA::normalize(Xt)
  Xt <- XN$Xc
  # Augment the data
  augX <- matrix(0,0,dim(Xt)[2]+K-1)
  augY <- c()
  for(i in 1:(K-1)){
    qq <- i
    # Class 1
    inds <- max(1,qq-s+1):qq
    X1 <- matrix(0,0,dim(Xt)[2]+K-1)
    Y1 <- c()
    for(j in inds){
      inds1 <- which(Yt == j)
      X1 <- rbind(X1,cbind(Xt[inds1,],matrix(0,nrow = length(inds1),ncol = K-1)))
      Y1 <- c(Y1,rep(1,length(inds1)))
    }
    # Class  2
    inds <- (qq+1):min(K,qq+s)
    X2 <- matrix(0,0,dim(Xt)[2]+K-1)
    Y2 <- c()
    for(j in inds){
      inds2 <- which(Yt == j)
      X2 <- rbind(X2,cbind(Xt[inds2,],matrix(0,nrow = length(inds2),ncol = K-1)))
      Y2 <- c(Y2,rep(2,length(inds2)))
    }
    X <- rbind(X1,X2)
    X[,dim(Xt)[2]+qq] <- rep(h,dim(X)[1]) # TRY TO CHANGE TO not qq-1
    Y <- c(Y1,Y2)
    augX <- rbind(augX,X)
    augY <- c(augY,Y)
  }
  augY <- factor(augY)
  # Train the model
  res <- accSDA::ASDA(Xt = augX, Yt = augY, Om= Om, control=list(ordinal=TRUE),...)
  res$XN <- XN
  class(res) <- 'ordASDA'
  res$call <- match.call()
  res$K <- K
  res$h <- h
  # Change the output...
  return(res)
}

#' @export
ordASDA.data.frame <- function(Xt, ...){
  res <- ordASDA.default(Xt = structure(data.matrix(Xt), class = "matrix"), ...)
  cl <- match.call()
  cl[[1L]] <- as.name("ordASDA")
  res$call <- cl
  res
}

#' @export
ordASDA.matrix <- function(Xt, ...){
  res <- ordASDA.default(Xt, ...)
  cl <- match.call()
  cl[[1L]] <- as.name("ASDA")
  res$call <- cl
  res
}

#' Predict method for ordinal sparse discriminant analysis
#'
#' Predicted values based on fit from the function \code{\link{ordASDA}}. This
#' function is used to classify new observations based on their explanatory variables/features.
#' There is no need to normalize the data, the data is normalized based on the normalization
#' data from the ordASDA object.
#'
#' @param object Object of class ordASDA. This object is returned from the function \code{\link{ordASDA}}.
#' @param newdata A matrix of new observations to classify.
#' @param ... Arguments passed to \code{\link[MASS]{predict.lda}}.
#' @return A vector of predictions.
#' @seealso \code{\link{ordASDA}}
#' @examples
#'     set.seed(123)
#'
#'     # You can play around with these values to generate some 2D data to test one
#'     numClasses <- 15
#'     sigma <- matrix(c(1,-0.2,-0.2,1),2,2)
#'     mu <- c(0,0)
#'     numObsPerClass <- 5
#'
#'     # Generate the data, can access with train$X and train$Y
#'     train <- accSDA::genDat(numClasses,numObsPerClass,mu,sigma)
#'     test <- accSDA::genDat(numClasses,numObsPerClass*2,mu,sigma)
#'
#'     # Visualize it, only using the first variable gives very good separation
#'     plot(train$X[,1],train$X[,2],col = factor(train$Y),asp=1,main="Training Data")
#'
#'     # Train the ordinal based model
#'     res <- accSDA::ordASDA(train$X,train$Y,s=2,h=1, gam=1e-6, lam=1e-3)
#'     vals <- predict(object = res,newdata = test$X) # Takes a while to run ~ 10 seconds
#'     sum(vals==test$Y)/length(vals) # Get accuracy on test set
#'     plot(test$X[,1],test$X[,2],col = factor(test$Y),asp=1,
#'          main="Test Data with correct labels")
#'     plot(test$X[,1],test$X[,2],col = factor(vals),asp=1,
#'          main="Test Data with predictions from ordinal classifier")
#' @rdname predict.ordASDA
#' @export
predict.ordASDA <- function(object, newdata = NULL, ...){
  newdata <- accSDA::normalizetest(newdata,object$XN)
  K <- object$K
  h <- object$h
  pred <- c()
  class(object) <- "ASDA"
  # Now we make all possible copies of the data, predict for
  # each and count the number of 2's
  for(i in 1:dim(newdata)[1]){
    classes <- c()
    for(j in 1:(K-1)){
      extra <- rep(0,K-1)
      extra[j] <- h
      obs <- matrix(c(newdata[i,],extra),nrow = 1)
      classes <- c(classes,stats::predict(object,newdata = obs)$class)
    }
    #print(i)
    #print(classes)
    k <- sum(classes=="2")+1
    pred <- c(pred,k)
  }
  return(pred)
}

#' Generate data for ordinal examples in the package
#'
#' Given the parameters, the function creates a dataset for testing the ordinal functionality
#' of the package. The data is samples from multivariate Gaussians with different means, where
#' the mean varies along a sinusoidal curve w.r.t. the class label.
#'
#' @param numClasses Positive integer specifying the number of classes for the
#'                   dataset.
#' @param numObsPerClass Number of observations sampled per class.
#' @param mu Mean of the first class.
#' @param sigma 2 by 2 covariance matrix
#'
#' @return \code{genDat} Returns a list with the following attributes:
#'     \describe{
#'       \item{X}{A matrix with two columns and \code{numObsPerClass}*\code{numClasses} rows.}
#'       \item{Y}{Labels for the rows of \code{X}.}
#'       }
#' @author Gudmundur Einarsson
#' @details
#' This function is used to demonstrate the usage of the ordinal classifier.
#' @seealso \code{\link{ordASDA}}
#' @examples
#' set.seed(123)
#'
#'     # You can play around with these values to generate some 2D data to test one
#'     numClasses <- 15
#'     sigma <- matrix(c(1,-0.2,-0.2,1),2,2)
#'     mu <- c(0,0)
#'     numObsPerClass <- 5
#'
#'     # Generate the data, can access with train$X and train$Y
#'     train <- accSDA::genDat(numClasses,numObsPerClass,mu,sigma)
#'     test <- accSDA::genDat(numClasses,numObsPerClass*2,mu,sigma)
#'
#'     # Visualize it, only using the first variable gives very good separation
#'     plot(train$X[,1],train$X[,2],col = factor(train$Y),asp=1,main="Training Data")
#'
#' @rdname genDat
#' @export genDat
genDat <- function(numClasses,numObsPerClass,mu,sigma){
  Xt <- matrix(0,0,2)
  Yt <- c()
  for(i in 1:numClasses){
    Xt <- rbind(Xt,MASS::mvrnorm(numObsPerClass,mu = mu,Sigma = sigma))
    Yt <- c(Yt,rep(i,numObsPerClass))
    mu[1] <- mu[1]+5
    mu[2] <- 20*sin((i/numClasses)*4*pi)
    #mu[2] <- 20*rnorm(1)
  }
  return(list(X=Xt,Y=Yt))
}
