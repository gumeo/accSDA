#' Cross-validation of sparse zero variance discriminant analysis
#'
#' Applies alternating direction methods of multipliers to solve sparse
#' zero variance discriminant analysis.
#'
#' @param X n by p data matrix, variables should be scaled to by sd
#' @param Y n by K indicator matrix.
#' @param folds number of folds to use in K-fold cross-validation.
#' @param gams Number of regularly spaced regularization parameters to try in [0,1]*max_gamma.
#'        See details for how max_gamma is computed in the function.
#' @param beta Augmented Lagrangian parameter. Must be greater than zero.
#' @param D Penalty dictionary basis matrix.
#' @param q Desired number of discriminant vectors.
#' @param maxits Number of iterations to run ADMM algorithm.
#' @param tol Stopping tolerances for ADMM, must have tol$rel and tol$abs.
#' @param ztol Rounding tolerance for truncating entries to 0.
#' @param feat Maximum fraction of nonzero features desired in validation scheme.
#' @param penalty Controls whether to apply reweighting of l1-penalty (using sigma = within-class std devs)
#' @param quiet toggles between displaying intermediate statistics.
#' @param ... Parameters passed to SZVD.default.
#' @return \code{SZVDcv} returns an object of \code{\link{class}} "\code{SZVDcv}"
#'        including a list with the named components \code{DVs} and \code{gambest}.
#'        Where \code{DVs} are the discriminant vectors for the best l1 regularization
#'        parameter and \code{gambest} is the best regularization parameter found
#'        in the cross-validation.
#' @seealso Used by: \code{\link{SZVDcv}}.
#' @examples
#'   P <- 150 # Number of variables
#'   N <- 20 # Number of samples per class
#'
#'   # Mean for classes, they are zero everywhere except the first 3 coordinates
#'   m1 <- rep(0,P)
#'   m1[1] <- 3
#'
#'   m2 <- rep(0,P)
#'   m2[2] <- 3
#'
#'   m3 <- rep(0,P)
#'   m3[3] <- 3
#'
#'   # Sample dummy data
#'   Xtrain <- rbind(MASS::mvrnorm(n=N,mu = m1, Sigma = diag(P)),
#'                  MASS::mvrnorm(n=N,mu = m2, Sigma = diag(P)),
#'                 MASS::mvrnorm(n=N,mu = m3, Sigma = diag(P)))
#'
#'
#'   # Generate the labels
#'   Ytrain <- cbind(c(rep(1,N),rep(0,2*N)),
#'                   c(rep(0,N),rep(1,N),rep(0,N)),
#'                   c(rep(0,2*N),rep(1,N)))
#'
#'   # Normalize the data
#'   Xt <- accSDA::normalize(Xtrain)
#'   Xtrain <- Xt$Xc
#'
#'   # Train the classifier and increase the sparsity parameter from the default
#'   # so we penalize more for non-sparse solutions.
#'   res <- accSDA::SZVD_kFold_cv(Xtrain,Ytrain,folds=2,gams=2,beta=2.5,q=1, D=diag(P),
#'                               maxits=50,tol=list(abs=1e-2,rel=1e-2),
#'                               ztol=1e-3,feat=0.3,quiet=FALSE,penalty=TRUE)
#' @details
#' Add how max_gamma is calculated from the ZVD solution.
#' This function might require a wrapper similar to ASDA.
#' @rdname SZVD_kFold_cv
#' @export SZVD_kFold_cv
SZVD_kFold_cv <- function(X, ...) UseMethod("SZVD_kFold_cv", X)

#' @return \code{NULL}
#' @export
#' @rdname SZVD_kFold_cv
#' @method SZVD_kFold_cv default
SZVD_kFold_cv.default <- function(X, Y, folds, gams, beta, D, q, maxits, tol, ztol, feat, penalty, quiet, ...){

  # Try to put everything in the style of the old function to make the code for experiments more convenient.
  # Also, let's choose the best parameter in the same way as for the other methods, meaning that we need
  # a feats parameter, setting a threshold for the maximum fraction of nonzero values in the discriminant vectors.
  ## Initialization

  # Get dimensions of input matrices
  n <- dim(X)[1]
  p <- dim(X)[2]
  K <- dim(Y)[2]
  num_gammas <- gams

  # If n is not divisible by K, duplicate some records for the sake of
  # cross validation.
  pad <- 0
  if(n %% folds > 0){
    pad <- ceiling(n/folds)*folds - n

    # Add the duplicates, such that number of data points is
    # divisible by the number of folds
    X <- rbind(X,X[1:pad,])
    Y <- rbind(Y,Y[1:pad,])
  }

  # Get the new number of rows
  n <- dim(X)[1]

  # Randomly permute rows of X
  prm <- sample(1:n,n,replace=FALSE)
  X <- X[prm,]
  Y <- Y[prm,]

  # Make Atrain
  Atrain <- cbind(Y%*%matrix(1:K,K,1),X)

  ###
  # Initialization of cross-validation indices
  ###

  # Number of validation samples
  nv <- n/folds

  # Initial validation indices
  vinds <- 1:nv

  # Initial training indices
  tinds <- (nv+1):n

  # Validation scores
  scores <- q*p*matrix(1,nrow = folds, ncol = num_gammas)

  # Misclassification rate for each classifier
  mc <- matrix(0,nrow = folds, ncol = num_gammas)

  for(f in 1:folds){
    if(!quiet){
      print("-------------------------------------------")
      print(paste("Fold number:",f))
      print("-------------------------------------------")
    }
    ## Initialization

    # Extract X and Y data
    Xt <- X[tinds,]
    At <- Atrain[tinds,]

    # Extract validation data
    Av <- Atrain[vinds,]

    # Get dimensions of training matrices.
    nt <- dim(Xt)[1]
    p <- dim(Xt)[2]

    resCV <- SZVDcv(Atrain = At,
                    Aval = Av,
                    k=K,
                    num_gammas = num_gammas,
                    g_mults = c(0,1),
                    D = diag(p),
                    sparsity_pen = feat,
                    scaling = FALSE,
                    penalty = penalty,
                    beta = beta,
                    tol = tol,
                    ztol = ztol,
                    maxits = maxits,
                    quiet = quiet)

    scores[f,] <- resCV$scores
  } # End folds

  ###
  # Find the best solution
  ###

  # Average CV scores
  avg_score <- colMeans(scores)

  # Choose lambda with best average score
  gbest <- which.min(avg_score)

  ###
  # Solve with lambda = lam(lbest)
  ###
  print(paste("Finished Training: best gamma ind =", gbest))

  # Use the full training set to obtain parameters
  Xt <- X[1:(n-pad),]
  Yt <- Y[1:(n-pad),]
  Atrain <- cbind(Yt%*%matrix(1:K,K,1),Xt)

  ###
  # Generate the gammas
  ###

  # Call ZVD function to solve the unpenalized problem.
  w0 = ZVD(Atrain, scaling=FALSE, get_DVs=TRUE)

  # Extract scaling vector for weighted l1 penalty and diagonal penalty matrix.
  if (penalty==TRUE){ # scaling vector is the std deviations of each feature.
    s = sqrt(diag(w0$W))
  }  else if(penalty==FALSE){ # scaling vector is all-ones (i.e., no scaling)
    s = rep(1, times=p)
  }
  w0$s = s

  ##################################################################################
  ## Compute range of sensible parameter values.
  ##################################################################################

  ## Normalize B (divide by the spectral norm)
  if (dim(w0$B)[2]==1){
    w0$B = w0$B/norm(w0$B, type='f')
  }  else{
    w0$B = (w0$B + t(w0$B))/eigen((w0$B + t(w0$B)), symmetric=TRUE, only.values=TRUE)$values[1]
  }

  # Compute ratio of max gen eigenvalue and l1 norm of the first ZVD to get "bound" on gamma.
  if (dim(w0$B)[2]==1){
    max_gamma =  (t(w0$dvs)%*%w0$B)^2/sum(abs(s*(D %*% w0$dvs)))
  }else{
    max_gamma = apply(w0$dvs, 2, function(x){(t(x) %*% w0$B %*% x)/sum(abs(s*(D%*%x)))})
  }

  # Generate range of gammas to choose from.
  g_mults <- c(0,1) # hardcoded for now
  gammas = sapply(max_gamma, function(x){seq(from=g_mults[1]*x, to=g_mults[2]*x, length=num_gammas)})
  gambest = gammas[gbest,]
  ###

  # Loop until nontrivial solution is found
  trivsol <- TRUE
  while(trivsol){
    szvdObj <- SZVD(Atrain, gambest, D, penalty, FALSE, tol, maxits, beta,TRUE)
    DVs <- szvdObj$DVs

    # Round small entried to zero
    DVs <- DVs*(as.matrix((abs(DVs)-ztol)>0)*1)

    # Check for trivial solution or Infs
    if(any(!is.finite(DVs))){
      # If infinite values we check the next solution
      #gbest <- gbest - 1
      gbest <- gbest + 1
      if(gbest > dim(gammas)[1]){
        gbest <- gbest - 1
        break
      }
      gambest <- gammas[gbest,]
    }else if(sum(DVs != 0) == 0){
      # If trivial solution
      gbest <- gbest + 1
      if(gbest > dim(gammas)[1]){
        gbest <- gbest - 1
        break
      }
      gambest <- gammas[gbest,]
    }else{
      trivsol <- FALSE
    }
  }

  # Create an object of class SDAPcv to return, might add more to it later
  retOb <- structure(
    list(call = match.call(),
         DVs = DVs,
         gambest = gambest,
         bestInd = gbest),
    class = "SZVD_kFold_cv")

  return(retOb)
}

#' @export
SZVD_kFold_cv.matrix <- function(X, ...){
  res <- SZVD_kFold_cv.default(X, ...)
  cl <- match.call()
  cl[[1L]] <- as.name("SZVD_kFold_cv")
  res$call <- cl
  res
}
