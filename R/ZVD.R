#' Zero Variance Discriminant Analysis
#'
#' Implements the ZVD algorithm to solve dicriminant vectors.
#'
#' @param A Matrix, where first column corresponds to class labels.
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
#'   \item{\code{class_obs}}{Matrices of observations of each class.}
#' }
#' @seealso Used by: \code{\link{SZVDcv}}.
#' @details
#' This function should potentially be made internal for the release.
#' @rdname ZVD
#' @export ZVD
ZVD <- function (A, ...) UseMethod("ZVD",A)

#' @return \code{NULL}
#'
#' @rdname ZVD
#' @method ZVD default
ZVD.default <- function(A, scaling = FALSE, get_DVs = FALSE){
  classes = factor(A[,1])
  X = as.matrix(data.frame(A[,2:dim(A)[2]]))

  # Get input dimensions.
  n = dim(X)[1]
  p = dim(X)[2]

  # Center the data.
  mu = colMeans(X)
  X = X - (rep(1, times= n )) %*% t(mu)

  ## Scale the data so each feature has variance equal to 1.
  if (scaling){
    # Compute standard deviation of each feature.
    sig = apply(X=X,MARGIN=2, FUN = function(y){
      if(sd(y)==0) {return(1)}
      else{return(sd(y))}
    }
    )

    # Divide each feature by its standard deviation to ensure each feature has variance equal to one.
    X = X%*% diag(1/sig)
  }





  #########################################################
  #### Extract classes from the observations.
  ########################################################

  # Initialize observations as an empty list.
  class_obs = list()

  # Extract class labels.
  labels = levels(classes)

  # Get number of classes.
  K = length(labels)

  # Initialize matrix of within-class means.
  classMeans = matrix(0, nrow = p, ncol = K)

  # Initialize within-class covariance.
  W = matrix(0, nrow = p, ncol = p)

  # Initialize between-class covariance (if more than 2 classes)
  if (K>2){
    B = matrix(0, nrow = p, ncol = p)
  }

  # DIAG: sizes.
  sizes = rep(0, K)
  #######################################################
  # For each class, make an object in the list containing only the observations in that class
  # and update the between and within-class sample covariances.
  #######################################################
  for (i in 1: K){

    # Find indices of observations in the current class.
    class_inds = (classes == labels[i])

    # Make new object in class_obs corresponding to the observations in this class.
    class_obs[[i]] = X[(class_inds==T) , ]

    # Get number of observations in this class.
    ni = dim(class_obs[[i]])[1]
    sizes[i] = ni

    # Compute within-class means.
    classMeans[,i] = colMeans(class_obs[[i]])

    # Compute within-class covariance (from the formula)
    for (j in 1:ni){
      xj = as.matrix(class_obs[[i]][j,])

      # Update W.
      W = W + (xj - classMeans[,i]) %*% t(xj - classMeans[,i])

    }

    # Update B (K>2 case)
    if (K>2){
      B = B + ni*classMeans[,i] %*% t(classMeans[,i])
    }

  }

  # Symmetrize W.
  W = (W + t(W))/2

  # Compute B (K=2 case) (Store as a p by 1 matrix/vector)
  if (K==2){
    B = as.matrix((classMeans[,1]-classMeans[,2]))
  }



  #######################################################
  ## Find ZVDs (if wanted)
  #######################################################
  if (get_DVs){
    # Find the null-vectors of W.
    S = eigen(W, symmetric=TRUE,)
    ds = sort(S$values, index.return=TRUE, decreasing=TRUE)
    V = as.matrix(S$vectors[,ds$ix])
    zeros = (ds$x < 1e-6)
    N = as.matrix(V[,zeros])

    # Compute the Zero Variance Discriminants.
    if (K==2){
      # Compute max eigenvector of N'*B*N
      w = N%*%t(N)%*%(classMeans[,1] - classMeans[,2])
      w = w/norm(as.matrix(w), 'F')
    }else{
      # Compute K-1 nontrivial eigenvectors of N'*B*N.
      w = svd(t(N) %*% B %*%N, nu=(K-1), nv= (K-1))$v[,1:K-1]
      # Project back to the original space.
      w = N %*% w
    }


  }

  #######################################################
  # Prep output.
  #######################################################

  # Output scaling/centring terms (if used).
  if (scaling){
    mu = list(mu=mu, sig=sig)
  }

  # Output list.
  if (get_DVs){
    ZVD = list(dvs=w, B=B, W=W, N=N, mu=mu, means=classMeans, k=K, labels=classes, obs=X, class_obs=class_obs, sizes=sizes)
  } else{
    ZVD = list(B=B, W=W, mu=mu, means=classMeans, k=K, labels=classes, obs=X, class_obs=class_obs, sizes=sizes)
  }

  ## Return output.
  return(ZVD)
}

#' @export
ZVD.matrix <- function(A, ...){
  res <- ZVD.default(A, ...)
  res
}
