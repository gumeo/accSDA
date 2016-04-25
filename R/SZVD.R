#' Sparse Zero Variance Discriminant Analysis
#'
#' Applies SZVD heuristic for sparse zero-variance discriminant 
#' analysis to given training set.
#'
#' @param train Data matrix where first column is the response class.
#' @param gamma Regularization paramter controlling l1-penalty.
#' @param D Penalty dictionary basis matrix.
#' @param penalty Controls wether to apply reweighting of l1-penalty 
#'        (using sigma = within-class std devs)
#' @param scaling Logical indicating whether to scale data such that each
#'        feature has variance 1.
#' @param tol Stopping tolerances for ADMM algorithm, 
#'        must include tol$rel and tol$abs.
#' @param maxits Maximum number of iterations used in the ADMM algorithm.
#' @param beta Parameter for augmented Lagrangian term in the ADMM algorithm.
#' @return \code{SZVD} returns an object of \code{\link{class}} "\code{SZVD}" including a list
#' with the following named components
#'
#' \describe{
#'   \item{\code{DVs}}{Discriminant vectors.}
#'   \item{\code{its}}{Number of iterations required to find DVs.}
#'   \item{\code{pen_scal}}{Weights used in reweighted l1-penalty.}
#'   \item{\code{N}}{Basis for the null-space of the sample within-class covariance.}
#'   \item{\code{w0}}{unpenalized zero-variance discriminants (initial solutions) plus B and W, etc.}
#' }
#' @seealso Used by: \code{\link{SZVDcv}}.
#' @details
#' This function will currently solve as a standalone function in accSDA for time comparison.
#' A wrapper function like ASDA will be created to use the functionality of plots and such.
#' Maybe call it ASZDA. For that purpose the individual ZVD function will need to be implemented.
#' @rdname SZVD
#' @export SZVD
SZVD <- function (x, ...) UseMethod("SZVD")

#' @return \code{NULL}
#'
#' @rdname SZVD
#' @method SZVD default
SZVD.default <- function(train, gamma, D, penalty, scaling, tol, maxits, beta, quiet){
  #==============================================================================================
  # Preprocess the training set.
  #==============================================================================================
    
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Get covariance matrices and initial solutions.
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # Get dimensions of the training set
  p <- dim(train)[2]
  p <- p-1
  
  # Call ZVD to process the training data
  get_DVs <- TRUE
  w0 <- ZVD(train, scaling, get_DVs)
  
  # Normalize B (divide by the spectral norm)
  Btype <- dim(w0$B)[2]
  if(Btype == 1){
    w0$B <- w0$B/norm(w0$B, type = "2")
  } else{
    w0$B <- 0.5*(w0$B + t(w0$B))
    w0$B <- w0$B/norm(w0$B, type = "2")
  }
  
  # Extract scaling vector for weighted l1 penalty and diagonal penalty matrix
  if(penalty){
    # scaling vector is the standard deviations of each feature.
    s <- sqrt(diag(w0$W))
  } else{
    # scaling vector is all-ones (i.e., no scaling)
    s <- matrix(1,p,1)
  }
  
  #==============================================================================================
  # Initialization for the algorithm.
  #==============================================================================================
  
  # Initialize output
  DVs <- matrix(0,p,w0$k-1)
  its <- matrix(0,w0$k-1,1)
  
  # Initialize objective matrix
  if(Btype == 1){
    B0 <- w0$B
  } else{
    B0 <- crossprod(w0$N,w0$B)%*%w0$N
    B0 <- 0.5*(B0 + t(B0))
  }
  
  # Initialize N
  N <- w0$N
  
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
    admmObj <- SZVD_ADMM(B0,  N, D, sols0, s,  gamma[i], beta, tol, maxits, quiet)
    x <- admmObj$x
    numIts <- admmObj$its
    
    # Save output
    DVs[,i] <- D%*%(N%*%x)
    its[i] <- numIts
    
    if(!quiet){
      print(paste("Found SZVD ", i, " after ", its[i], " iterations" ,sep=""))
    }
    
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Update N using QR.
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(i < (w0$k-1)){
      # Project columns of N onto orthogonal complement of Nx.
      x <- DVs[,i]
      x <- x/norm(x,type="2")
      
      Ntmp <- N - x%*%crossprod(x,N)
      
      # Call QR factorization to extract ortonormal basis for span(Ntmp)
      qrObj <- qr(Ntmp)
      Q <- qr.Q(qrObj)
      R <- qr.R(qrObj)
      
      # Extract nonzero rows of R to get columns of Q to use as new N
      R_rows <- (abs(diag(R)) > 1e-6)
      
      # use nontrivial columns of Q as updated N
      N <- Q[,R_rows]
      
      # Update B0 according to the new basis N
      B0 <- crossprod(N,w0$B)%*%N
      B0 <- 0.5*(B0 + t(B0))
    }
  }
  
  #==============================================================================================
  # Prep output.
  #==============================================================================================
  outOb <- structure(list(call = match.call(),
                          DVs = DVs,
                          its = its,
                          pen_scal = s,
                          N = N,
                          w0 = w0),
                     class = "SZVD")
  return(outOb)
}