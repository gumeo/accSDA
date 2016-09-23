#' Sparse Zero Variance Discriminant Analysis
#'
#' Applies SZVD heuristic for sparse zero-variance discriminant
#' analysis to given training set.
#'
#' @param train Data matrix where first column is the response class.
#' @param gamma Set of regularization parameters controlling l1-penalty.
#' @param D dictionary/basis matrix.
#' @param penalty Controls whether to apply reweighting of l1-penalty (using sigma = within-class std devs)
#' @param scaling Logical indicating whether to scale data such that each
#'        feature has variance 1.
#' @param tol Stopping tolerances for ADMM algorithm,
#'        must include tol$rel and tol$abs.
#' @param maxits Maximum number of iterations used in the ADMM algorithm.
#' @param penalty Controls whether to apply reweighting of l1-penalty (using sigma = within-class std devs).
#' @param quiet Print intermediate outpur or not.
#' @return \code{SZVD} returns an object of \code{\link{class}} "\code{SZVD}" including a list
#' with the following named components:
#'
#' \describe{
#'   \item{\code{DVs}}{Discriminant vectors.}
#'   \item{\code{its}}{Number of iterations required to find DVs.}
#'   \item{\code{pen_scal}}{Weights used in reweighted l1-penalty.}
#'   \item{\code{N}}{Basis for the null-space of the sample within-class covariance.}
#'   \item{\code{means}}{Training class-means.}
#'   \item{\code{mus}}{Training meand and variance scaling/centring terms.}
#'   \item{\code{w0}}{unpenalized zero-variance discriminants (initial solutions) plus B and W, etc.}
#' }
#' @seealso Used by: \code{\link{SZVDcv}}.
#' @details
#' This function will currently solve as a standalone function in accSDA for time comparison.
#' A wrapper function like ASDA will be created to use the functionality of plots and such.
#' Maybe call it ASZDA. For that purpose the individual ZVD function will need to be implemented.
#' @rdname SZVD
#' @export SZVD
SZVD <- function (train, ...) UseMethod("SZVD",train)

#' @return \code{NULL}
#'
#' @rdname SZVD
#' @method SZVD default
SZVD.default <- function(train, gamma, D, penalty=TRUE, scaling=TRUE, tol = list(abs=1e-4, rel=1e-4), maxits = 2000, beta=1, quiet=TRUE){
  #==============================================================================================
  # Preprocess the training set.
  #==============================================================================================

  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Get covariance matrices and initial solutions.
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  # Get dimensions of the training set
  n <- dim(train)[1]
  p <- dim(train)[2]-1

  # Call ZVD to process the training data
  w0 <- ZVD(train, scaling = scaling, get_DVs = TRUE)

  # Normalize B (divide by the spectral norm)
  if(dim(w0$B)[2] == 1){
    w0$B <- w0$B/norm(w0$B, type = "f")
  } else{
    w0$B = 0.5*(w0$B + t(w0$B))/eigen(0.5*(w0$B + t(w0$B)), symmetric=TRUE, only.values=TRUE)$values[1]
  }

  # Force N to be a matrix
  w0$N <- as.matrix(w0$N)

  # Extract scaling vector for weighted l1 penalty and diagonal penalty matrix
  if(penalty == TRUE){
    # scaling vector is the standard deviations of each feature.
    s <- sqrt(diag(w0$W))
  } else{
    # scaling vector is all-ones (i.e., no scaling)
    s = rep(1, times=p)
  }

  #==============================================================================================
  # Initialize missing arguments
  #==============================================================================================

  # If gamma is missing, use ratio of objectives for the ZVDs to get "good" guess for gamma.
  if (missing(gamma)){
    if (dim(w0$B)[2]==1){
      gamma = (t(w0$B)%*%w0$dvs)^2/sum(abs(s*w0$dvs))
    }   else{
      gamma = apply(w0$dvs, 2, function(x){(t(x) %*% w0$B %*% x)/sum(abs(s*x))})
    }
    ## Scale gamma.
    gamma=0.75*gamma
  }

  # If dictionary D missing, use the identity matrix.
  if (missing(D)){
    D = diag(p)
  }

  #==============================================================================================
  # Initialization for the algorithm.
  #==============================================================================================

  # Initialize output.
  DVs = matrix(0, nrow = p, ncol = (w0$k-1))
  its = rep(0, times=(w0$k-1))

  # Initalize objective matrix
  if (dim(w0$B)[2]==1){
    B0 = w0$B
  }    else{
    B0 = t(w0$N) %*% w0$B %*% w0$N
    B0 = 0.5*(B0+t(B0))
  }

  # Initialize nullspace matrix.
  N = w0$N

  #==============================================================================================
  # Find the DVs sequentially using ADMM.
  #==============================================================================================

  for(i in 1:(w0$k-1)){
    # Initial solution
    sols0 <- list(x = t(N)%*%crossprod(D,w0$dvs[,i]),
                  y = w0$dvs[,i],
                  z = matrix(0,p,1))

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Call ADMM solver.
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    tmp = SZVD_ADMM(B = B0,  N = N, D=D, pen_scal=s,
                    gamma=gamma[i], beta=beta, tol=tol, maxits=maxits, quiet=TRUE)

    # Save output.
    DVs[,i] = matrix(D%*%N%*%tmp$x, nrow=p, ncol=1)
    its[i] = tmp$its

    if (quiet == FALSE){
      print(sprintf('Found SZVD # %g after %g its', i, its[i]) ,quote=FALSE)

    }

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Update N using QR.
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (i < (w0$k-1)) {
      # Project columns of N onto orthogonal complement of Nx.
      x = DVs[,i]
      x = x/norm(as.matrix(x), 'f')

      Ntmp = N - x %*% (t(x)%*%N)

      # Call QR factorization to extract orthonormal basis for span(Ntmp)
      QR = qr(Ntmp)

      # Extract nonzero rows of R to get columns of Q to use as new N.
      R_rows = (abs(diag(qr.R(QR))) > 1e-6)

      # Use nontrivial columns of Q as updated N.
      N = qr.Q(QR)[, R_rows]

      # Update B0 according to the new basis N.
      B0 = t(N) %*% w0$B %*% N
      B0 = 0.5*(B0+t(B0))
    } # End if.
  } # End for.

  #==============================================================================================
  # Prep output.
  #==============================================================================================
  return(structure(list(call = match.call(),
                        DVs = DVs,
                        its = its,
                        pen_scal = s,
                        N = N,
                        means = w0$means,
                        mus = w0$mu,
                        w0 = w0),
                   class = "SZVD"))
}

#' @export
SZVD.matrix <- function(train, ...){
  res <- SZVD.default(train, ...)
  res
}
