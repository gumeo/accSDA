#' Sparse Discriminant Analysis solved via Proximal Gradient
#'
#' Applies proximal gradient algorithm to
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
#' @param q Desired number of discriminant vectors.
#' @param PGsteps Maximum number if inner proximal gradient algorithm for finding beta.
#' @param PGtol Stopping tolerance for inner APG method.
#' @param maxits Number of iterations to run
#' @param tol Stopping tolerance for proximal gradient algorithm.
#' @param initTheta Initial first theta, default value is a vector of ones.
#' @param bt Boolean to indicate whether backtracking should be used, default false.
#' @param L Initial estimate for Lipshitz constant used for backtracking.
#' @param eta Scalar for Lipshitz constant.
#' @return \code{SDAP} returns an object of \code{\link{class}} "\code{SDAP}" including a list
#' with the following named components: (More will be added later to handle the predict function)
#'
#' \describe{
#'   \item{\code{call}}{The matched call.}
#'   \item{\code{B}}{p by q matrix of discriminant vectors.}
#'   \item{\code{Q}}{K by q matrix of scoring vectors.}
#' }
#' @seealso \code{SDAPcv}, \code{\link{SDAAP}} and \code{\link{SDAD}}
#' @keywords internal
SDAP <- function (x, ...) UseMethod("SDAP")

#' @return \code{NULL}
#'
#' @rdname SDAP
#' @method SDAP default
SDAP.default <- function(Xt, Yt, Om, gam, lam, q, PGsteps, PGtol, maxits, tol, initTheta, bt = FALSE, L, eta){

  # Read training data size
  n <- dim(Xt)[1]
  p <- dim(Xt)[2]
  K <- dim(Yt)[2]

  # Precompute repeatedly used matrix products
  A <- (crossprod(Xt) + gam*Om) # Elastic net coef matrix
  alpha <- 1/norm(A, type="2") # Step length in PGA
  #L <- 1/alpha
  L <- norm(diag(diag(Om*gam)),'I')+norm(Xt,'F')^2
  D <- (1/n)*(crossprod(Yt))
  R <- chol(D)

  # Initialize B and Q
  Q <- matrix(1,K,q)
  B <- matrix(0,p,q)
  #-----------------------------------------------------------
  # Alternating direction method to update (theta, beta)
  #-----------------------------------------------------------
  for(j in 1:q){
    ###
    # Initialization
    ###

    # Compute Qj (K by j, first j-1 scoring vectors, all-ones last col)
    Qj <- Q[,1:j]

    # Precompute Mj = I - Qj*Qj'*D
    Mj <- function(u){
      return(u-Qj%*%(crossprod(Qj,D%*%u)))
    }

    # Initialize theta
    theta <- matrix(stats::runif(K),nrow=K,ncol=1)
    theta <- Mj(theta)
    if(j == 1 & !missing(initTheta)){
      theta=initTheta/10
    }
    theta <- theta/as.numeric(sqrt(crossprod(theta,D%*%theta)))

    # In case we want to initialize the theta

    # Initialize beta
    beta <- matrix(0,p,1)

    ###
    # Alternating direction method to update (theta,beta)
    ###
    for(its in 1:maxits){
      # Compute coefficient vector for elastic net step
      d <- 2*crossprod(Xt,Yt%*%(theta/n))

      # Update beta using proximal gradient step
      b_old <- beta
      if(bt == FALSE){
        beta <- prox_EN(A, d, beta, lam, alpha, PGsteps, PGtol)
        beta <- beta$x
      }else{
        beta <- prox_ENbt(A, d, beta, lam, L, eta, PGsteps, PGtol)
        beta <- beta$x
      }


      # Update theta using the projected solution
      if(norm(beta, type="2") > 1e-12){
        b <- crossprod(Yt,Xt%*%beta)
        #y <- solve(t(R),b)
        #z <- solve(R,y)
        y <- forwardsolve(t(R),b)
        z <- backsolve(R,y)
        tt <- Mj(z)
        t_old <- theta
        theta <- tt/as.numeric(sqrt(crossprod(tt,D%*%tt)))

        # Progress
        db <- norm(beta-b_old)/norm(beta, type="2")
        dt <- norm(theta-t_old)/norm(theta, type="2")
      } else{
        # Update b and theta
        beta <- beta*0
        theta <- theta*0
        db <- 0
        dt <- 0
      }
      # Check convergence
      if(max(db,dt)<tol){
        # Converged
        break
      }
    }
    # Make the first argument be positive, this is to make the results
    # more reproducible and consistent.
    if(theta[1] < 0){
      theta <- (-1)*theta
      beta <- (-1)*beta
    }

    # Update Q and B
    Q[,j] <- theta
    B[,j] <- beta
  }
  #Return B and Q in a SDAP object
  retOb <- structure(
    list(call = match.call(),
         B = B,
         Q = Q),
    class = "SDAP")
  return(retOb)
}
