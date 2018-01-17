#' Sparse Discriminant Analysis solved via Accelerated Proximal Gradient
#'
#' Applies accelerated proximal gradient algorithm to
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
#' @param selector Vector to choose which parameters in the discriminant vector will be used to calculate the
#'                 regularization terms. The size of the vector must be *p* the number of predictors. The
#'                 default value is a vector of all ones. This is currently only used for ordinal classification.
#' @param initTheta Option to set the initial theta vector, by default it is a vector of all ones
#'                  for the first theta.
#' @param bt Boolean to indicate whether backtracking should be used, default false.
#' @param L Initial estimate for Lipshitz constant used for backtracking.
#' @param eta Scalar for Lipshitz constant.
#' @return \code{SDAAP} returns an object of \code{\link{class}} "\code{SDAAP}" including a list
#' with the following named components: (More will be added later to handle the predict function)
#'
#' \describe{
#'   \item{\code{call}}{The matched call.}
#'   \item{\code{B}}{p by q matrix of discriminant vectors.}
#'   \item{\code{Q}}{K by q matrix of scoring vectors.}
#' }
#' @seealso \code{SDAAPcv}, \code{\link{SDAP}} and \code{\link{SDAD}}
#' @keywords internal
SDAAP <- function(x, ...) UseMethod("SDAAP")


#' @return \code{NULL}
#'
#' @rdname SDAAP
#' @method SDAAP default
SDAAP.default <- function(Xt, Yt, Om, gam, lam, q, PGsteps, PGtol, maxits, tol, selector = rep(1,dim(Xt)[2]), initTheta, bt=FALSE, L, eta){
  # TODO: Handle Yt as a factor and generate dummy matrix from it

  # Get training data size
  nt <- dim(Xt)[1] # num. samples
  p <- dim(Xt)[2]  # num. features
  K <- dim(Yt)[2]  # num. classes

  if(length(selector) != dim(Xt)[2]){
    stop('The length of selector must be the same as that of Xt')
  }

  # Structure to store elements for matrix A used
  # later on, i.e. precomputed values for speed.
  A <- structure(list(
    flag = NA,
    gom = NA,
    X = NA,
    n = NA,
    A = NA
  ),
    class = "Amat"
  )

  # Check if Omega is diagonal
  if(norm(diag(diag(Om))-Om, type = "F") < 1e-15){
    A$flag <- 1
    A$gom <- gam*diag(Om)
    A$X <- Xt
    A$n <- nt
    A$A <- 2*(crossprod(Xt)/nt + gam*Om)
    alpha <- 1/(2*(norm(Xt, type="1")*norm(Xt, type="I")/nt + norm(diag(A$gom), type="I")))
  }else{
    A$flag <- 0
    A$A <- 2*(crossprod(Xt)/nt + gam*Om)
    alpha <- 1/(norm(A$A, type="F"))
  }

  D <- (1/nt)*crossprod(Yt)
  R <- chol(D)

  Q <- matrix(1,K,q)
  B <- matrix(0,p,q)

  #------------------------------------------------------
  # Alternating direction method to update (theta,beta)
  #------------------------------------------------------
  for(j in 1:q){
    ###
    # Initialization
    ###
    Qj <- Q[,1:j]

    Mj <- function(u){
      return(u-Qj%*%(crossprod(Qj,D%*%u)))
    }

    # Initialize theta
    theta <- matrix(stats::runif(K),nrow=K,ncol=1)
    theta <- Mj(theta)
    if(j == 1 & !missing(initTheta)){
      theta=initTheta
    }
    theta <- theta/as.numeric(sqrt(crossprod(theta,D%*%theta)))

    # Initialize beta
    beta <- matrix(0,p,1)

    for(its in 1:maxits){
      # Compute coefficient vector for elastic net step
      d <- 2*crossprod(Xt,Yt%*%(theta/nt))

      # Update beta using proximal gradient step
      b_old <- beta
      if(bt == FALSE){
        betaOb <- APG_EN2(A, d, beta, lam, alpha, PGsteps, PGtol, selector)
        beta <- betaOb$x
      }else{
        betaOb <- APG_EN2bt(A, d, beta, lam, L, eta, PGsteps, PGtol, selector)
        beta <- betaOb$x
      }

      # Update theta using the projected solution
      if(norm(beta, type="2") > 1e-12){
        # Update theta
        b <- crossprod(Yt,Xt%*%beta)
        #y <- solve(t(R),b)
        #z <- solve(R,y)
        y <- forwardsolve(t(R),b)
        z <- backsolve(R,y)
        tt <- Mj(z)
        t_old <- theta
        theta <- tt/sqrt(as.numeric(crossprod(tt, D)%*%tt))

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
  #Return B and Q in a SDAAP object
  retOb <- structure(
    list(call = match.call(),
         B = B,
         Q = Q),
    class = "SDAAP")
  return(retOb)
}

# TODO: Implement print, summary, predict and plot of SDAAP object
# Potentially also print for output from summary and predict.
# Skim over again to see if crossprod can be used for more speed.
