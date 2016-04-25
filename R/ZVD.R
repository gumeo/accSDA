#' Zero Variance Discriminant Analysis
#'
#' Implements the ZVD algorithm to solve dicriminant vectors.
#'
#' @param A Data matrix, where first column corresponds to class labels.
#' @param scaling Logical whether to rescale data so each feature has variance 1.
#' @param get_DVs Logical whether to obtain unpenalized zero-variance discriminant vectors.
#' @return \code{SZVDcv} returns an object of \code{\link{class}} "\code{ZVD}" 
#'        including a list with the following named components:
#' \describe{
#'   \item{\code{dvs}}{discriminant vectors (optional).}
#'   \item{\code{B}}{sample between-class covariance.}
#'   \item{\code{W}}{sample within-class covariance.}
#'   \item{\code{N}}{basis for the null space of the sample within-class covariance.}
#'   \item{\code{mu}}{training mean and variance scaling/centring terms}
#'   \item{\code{means}}{vectors of sample class-means.}
#'   \item{\code{k}}{number of classes in given data set.}
#'   \item{\code{labels}}{list of classes.}
#'   \item{\code{obs}}{matrix of data observations.}
#' }
#' @seealso Used by: \code{\link{SZVDcv}}.
#' @details
#' This function should potentially be made internal for the release.
#' @rdname ZVD
#' @export ZVD
ZVD <- function (A, ...) UseMethod("ZVD")

#' @return \code{NULL}
#'
#' @rdname ZVD
#' @method ZVD default
ZVD.default <- function(A, scaling, get_DVs){
  #
  # HERE WE NEED A DESCRIPTION
  # Use Roxygen2 to create the desired documentation
  #
  # TODO: Currently the labels are the first column of A
  #       this needs to be changed to be more general.
  # NOTE1: This function requires normalize from sparseLDA
  # Note2: Scaling is optional in the matlab version, maybe
  #        make a flag for normalization, i.e. center and scaling
  # Note3: Classes are assumed to be numeric values from 1 to K
  #        We need to be able to deal with factors in data frames
  #        as well.
  # Note4: Requires the function nullSp define in the same folder
  # Note5: eigs function from rARPACK used

  ###
  # Preprocessing
  ###

  # Extract class labels
  classes <- A[,1]

  # Extract observations
  X <- A[,2:(dim(A)[2])]

  # Normalize the training data
  Xnorm <- sparseLDA::normalize(X)
  X <- Xnorm$Xc
  n <- dim(X)[1]
  p <- dim(X)[2]

  # Get number of classes
  K <- length(unique(as.factor(classes)))

  # Initialize matrix of within-class means
  classMeans <- matrix(0, p, K)

  # Ininitialize within-class covariance
  W <- matrix(0, p, p)

  # Initialize between-class covariance (if more than 2 classes)
  if(K > 2){
    B <- matrix(0, p, p)
  }

  class_sizes <- matrix(0, K, 1)

  #-----------------------------------------------------------------------
  # Compute within-class scatter matrix
  #-----------------------------------------------------------------------
  for(i in 1:K){
    # Find indices of observations in the current class
    class_inds <- (classes == i)

    # Get the objects of this class
    class_obs <- X[class_inds,]

    # Get number of observations in this class
    ni <- dim(class_obs)[1]
    class_sizes[i,1] <- ni

    # Compute within-class means
    classMeans[,i] <- colMeans(class_obs)

    # Compute sample class-scatter matrix
    Xc <- class_obs - matrix(1,ni,1)%*%t(classMeans[,i])

    # Update W
    W <- W + t(Xc)%*%Xc
  }

  # Symmetrize W
  W <- (W + t(W))/2

  #-----------------------------------------------------------------------
  # Compute B
  #-----------------------------------------------------------------------
  if(K == 2){
    # (K=2 case) (Store as a p by 1 matrix/vector)
    B <- classMeans[,1,drop = FALSE] - classMeans[,2,drop = FALSE]
  }else{
    # K > 2 case.

    # Make partition matrix
    Y <- matrix(0, n, K)
    for(i in 1:n){
      # Note that here we are doing dummy coding of the
      # factor variable, this can probably be done more
      # efficiently
      Y[i,classes[i]] <- 1
    }

    # Set scatter matrix B
    XY <- t(X)%*%Y
    B <- XY%*%(solve(t(Y)%*%Y,t(XY)))

    # Symmetrize B
    B <- (B + t(B))/2
  }

  #-----------------------------------------------------------------------
  # Find the null-vectors of W
  #-----------------------------------------------------------------------
  N <- nullSp(W) # The null space calculated as in matlab

  #-----------------------------------------------------------------------
  # Find ZVDs (if GetDVs == TRUE)
  #-----------------------------------------------------------------------
  if(get_DVs){
    if(K==2){
      # Compute max eigenvector of t(N)%*%B%*%N
      w <- N%*%(t(N)%*%B)
      w <- w/norm(w, type = "2")
    } else{
      # Compute K-1 nontrivial eigenvectors of t(N)%*%B%*%N
      # This depends on the rARPACK package
      # This is for speed, so we do not need to
      # calculate all eigenvectors, that is slow
      val <- rARPACK::eigs(t(N)%*%B%*%N, k = K-1, which = 'LM')
      w <- val$vectors

      # Project back to original space
      w <- N%*%w
    }
  }

  #-----------------------------------------------------------------------
  # Prepare output
  #-----------------------------------------------------------------------
  outOb <- structure(list(call = match.call(),
                          mu = Xnorm$mx,
                          sig = Xnorm$vx,
                          dvs = NA,
                          B = B,
                          W = W,
                          N = N,
                          means = classMeans,
                          k = K,
                          labels = classes,
                          obs = X,
                          ni = class_sizes),
                     class = "ZVD")
  if(get_DVs){
    outOb$dvs <- w
  }
  return(outOb)
}

#' @export
ZVD.matrix <- function(A, ...){
  res <- ZVD.default(A, ...)
  res
}
