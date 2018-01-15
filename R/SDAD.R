#' Sparse Discriminant Analysis solved via ADMM
#'
#' Applies alternating direction methods of multipliers algorithm to
#' the optimal scoring formulation of sparse discriminant analysis proposed
#' by Clemmensen et al. 2011.
#'
#' @param Xt n by p data matrix, (not a data frame, but a matrix)
#' @param Yt n by K matrix of indicator variables (Yij = 1 if i in class j).
#'     This will later be changed to handle factor variables as well.
#'     Each observation belongs in a single class, so for a given row/observation,
#'     only one element is 1 and the rest is 0.
#' @param Om p by p parameter matrix Omega in generalized elastic net penalty.
#' @param gam Regularization parameter for elastic net penalty.
#' @param lam Regularization parameter for l1 penalty, must be greater than zero.
#' @param mu Penalty parameter for augmented Lagrangian term, must be greater than zero.
#' @param q Desired number of discriminant vectors.
#' @param PGsteps Maximum number if inner proximal gradient algorithm for finding beta.
#' @param PGtol Two stopping tolerances for inner ADMM method, first is absolute tolerance, second is relative.
#' @param maxits Number of iterations to run
#' @param tol Stopping tolerance for proximal gradient algorithm.
#' @param selector Vector to choose which parameters in the discriminant vector will be used to calculate the
#'                 regularization terms. The size of the vector must be *p* the number of predictors. The
#'                 default value is a vector of all ones. This is currently only used for ordinal classification.
#' @param initTheta Initial first theta, default value is a vector of ones.
#' @return \code{SDAD} returns an object of \code{\link{class}} "\code{SDAD}" including a list
#' with the following named components: (More will be added later to handle the predict function)
#' \describe{
#'   \item{\code{call}}{The matched call.}
#'   \item{\code{B}}{p by q matrix of discriminant vectors.}
#'   \item{\code{Q}}{K by q matrix of scoring vectors.}
#' }
#' @seealso \code{SDADcv}, \code{\link{SDAAP}} and \code{\link{SDAP}}
#' @keywords internal
SDAD <- function (x, ...) UseMethod("SDAD")

#' @return \code{NULL}
#'
#' @rdname SDAD
#' @method SDAD default
SDAD.default <- function(Xt, Yt, Om, gam, lam, mu, q, PGsteps, PGtol, maxits, tol, selector = rep(1,dim(Xt)[2]), initTheta){
  # TODO: Handle Yt as a factor and generate dummy matrix from it

  ###
  # Initialize training sets etc
  ###
  # Get training data size
  nt <- dim(Xt)[1] # num. samples
  p <- dim(Xt)[2]  # num. features
  K <- dim(Yt)[2]  # num. classes

  # Check if Om is diagonal. If so, use matrix inversion lemma in linear
  # system solves.
  if(norm(diag(diag(Om))-Om,type = "F") < 1e-15){
    # Flag to use Sherman-Morrison-Woodbury to translate to
    # smaller dimensional linear system solves.
    SMW <- 1

    # Easy to invert diagonal part of Elastic net coefficient matrix.
    M <- mu*diag(p) + 2*gam*Om
    Minv = 1/diag(M)

    # Cholesky factorization for smaller linear system.
    RS = chol(diag(nt) + 2*Xt%*%diag(Minv)%*%t(Xt)/nt);
  } else{ # Use Cholesky for solving linear systems in ADMM step
    # Flag to not use SMW
    SMW <- 0
    A <- mu*diag(p) + 2*(crossprod(Xt) + gam*Om)
    R2 <- chol(A)
  }

  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Matrices for theta update.
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  D <- (1/nt)*(crossprod(Yt))
  R <- chol(D) # Cholesky factorization of D.

  # Initialize B and Q.
  Q <- matrix(1,K,q)
  B <- matrix(0,p,q)

  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Call Alternating Direction Method to solve SDA.
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # For j=1,2,..., q compute the SDA pair (theta_j, beta_j).
  for(j in 1:q){
    ###
    # Initialization
    ###
    # Compute Qj (K by j, first j-1 scoring vectors, all-ones last col)
    Qj <- Q[,1:j]

    # Precompute Mj = I - Qj*Qj'*D.
    Mj <- function(u){
      return(u - Qj%*%(crossprod(Qj,D%*%u)))
    }

    # Initialize theta
    theta <- matrix(stats::runif(K),nrow=K,ncol=1)
    theta <- Mj(theta)
    if(j == 1 & !missing(initTheta)){
      theta=initTheta
    }
    theta <- theta/as.numeric(sqrt(crossprod(theta,D%*%theta)))

    # Initialize coefficient vector for elastic net step
    d <- 2*crossprod(Xt,Yt%*%theta)

    # Initialize beta
    if(SMW == 1){
      btmp <- Xt%*%(Minv*d)/nt
      #beta <- (Minv*d) - 2*Minv*(crossprod(Xt,solve(RS,solve(t(RS),btmp))))
      beta <- (Minv*d) - 2*Minv*(crossprod(Xt,backsolve(RS,forwardsolve(t(RS),btmp))))
    }else{
      #beta <- solve(R2,solve(t(R2),d))
      beta <- backsolve(R2,forwardsolve(t(R2),d))
    }

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Alternating direction method to update (theta, beta)
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++
    for(its in 1:maxits){
      # Update beta using alternating direction method of multipliers.
      b_old <- beta

      if(SMW == 1){
        # Use SMW-based ADMM
        betaOb <- ADMM_EN_SMW(Minv, Xt, RS, d, beta, lam, mu, PGsteps, PGtol, TRUE, selector)
        beta <- betaOb$y
      } else{
        betaOb <- ADMM_EN2(R2, d, beta, lam, mu, PGsteps, PGtol, TRUE, selector)
        beta <- betaOb$y
      }

      # Update theta using the projected solution
      if(norm(beta, type="2") > 1e-15){
        # Update theta
        b <- crossprod(Yt, Xt%*%beta)
        #y <- solve(t(R),b)
        #z <- solve(R,y)
        y <- forwardsolve(t(R),b)
        z <- backsolve(R,y)
        tt <- Mj(z)
        t_old <- theta
        theta <- tt/sqrt(as.numeric(crossprod(tt,D%*%tt)))

        # Update changes
        db <- norm(beta-b_old, type="2")/norm(beta, type="2")
        dt <- norm(theta-t_old, type="2")/norm(theta, type="2")
      } else{
        # Update b and theta
        beta <- beta*0
        theta <- theta*0
        db <- 0
        dt <- 0
      }

      if(max(db,dt) < tol){
        # Converged
        break
      }

    }
    # Update Q and B
    Q[,j] <- theta
    B[,j] <- beta
  }
  #Return B and Q in a SDAAP object
  retOb <- structure(
    list(call = match.call(),
         B = B,
         Q = Q),
    class = "SDAD")
  return(retOb)
}
