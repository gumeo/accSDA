#' Accelerated Proximal Gradient on l1 regularized quadratic program
#'
#' Applies accelerated proximal gradient algorithm to the l1-regularized quadratic program
#' (with rank reduced Omega inside A)
#' \deqn{f(\mathbf{x}) + g(\mathbf{x}) = \frac{1}{2}\mathbf{x}^TA\mathbf{x} - d^T\mathbf{x} + \lambda |\mathbf{x}|_1}{f(x) + g(x) = 0.5*x^T*A*x - d^T*x + lambda*|x|_l1}
#'
#' @param A Object containing everythign needed for calculating A,
#'        X and Omega are factored.
#' @param d nx1 dimensional column vector.
#' @param lam Regularization parameter for l1 penalty, must be greater than zero.
#' @param alpha Step length.
#' @param maxits Number of iterations to run
#' @param tol Stopping tolerance for proximal gradient algorithm.
#' @param selector Vector to choose which parameters in the discriminant vector will be used to calculate the
#'                 regularization terms. The size of the vector must be *p* the number of predictors. The
#'                 default value is a vector of all ones. This is currently only used for ordinal classification.
#' @return \code{APG_EN2} returns an object of \code{\link{class}} "\code{APG_EN2}" including a list
#' with the following named components
#'
#' \describe{
#'   \item{\code{call}}{The matched call.}
#'   \item{\code{x}}{Found solution.}
#'   \item{\code{k}}{Number of iterations used.}
#' }
#' @seealso Used by: \code{\link{SDAAP}} and the \code{SDAAPcv} cross-validation version.
#' @details
#' This function is used by other functions and should only be called explicitly for
#' debugging purposes.
#' @keywords internal
APG_EN2rr <- function(A, d, x0, lam, alpha,  maxits, tol, selector= rep(1,dim(x0)[1])){
  ###
  # Initialization
  ###

  x <- x0
  xold <- x

  # Get number of components of x
  p <- dim(x0)[1]

  # initial momentum coefficient
  t <- 1
  told <- 1

  # Objective function and gradient
  df <- function(x){
    2*(crossprod(A$X,A$X%*%x)+A$gom%*%t(A$gom)%*%x) - d
  }


  #-------------------------------------------------------------
  # Outer loop: Repeat until convergence or max # of iterations
  #-------------------------------------------------------------
  for(k in 0:maxits){
    # Compute gradient for differentiable part
    dfx <- df(x)

    #-------------------------------------------------------------
    # Compute disagreement between df and lam*sign(x) on supp(x)
    #-------------------------------------------------------------
    # Initialize error vector
    err <- matrix(0,p,1)
    # Initialize cardinality of support
    card <- 0

    # For each i, update error if i is in the support
    for(i in 1:p){
      if(abs(x[i]) > 1e-12){
        # Update cardinality
        card <- card + 1

        # Update error vector
        err[i] <- (-dfx[i]-lam*sign(x[i]))*(selector[i])
      }
    }

    #----------------------------------------------------------------
    # Check stopping criteria -df(x) in subdiff g(x).
    #   Need inf(df) < lam + tol, inf(err) < tol.
    #----------------------------------------------------------------
    if(max(norm(dfx,type="I")-lam,norm(err,type="I")) < tol*p){
      # CONVERGED!!!
      break
    } else{
      ###
      # Update x APG iteration
      ###
      # Compute extrapolation factor
      told <- t
      t <- (1 + sqrt(1+4*told^2))/2

      # Extrapolate using last two iterates
      y <- x+(told-1)/t*(x-xold)

      # Compute gradient at y
      dfy <- df(y)

      # Take proximal gradient step from y
      xold <- x
      xx <- sign(y-alpha*dfy)*pmax(abs(y-alpha*dfy) - lam*alpha*matrix(1,p,1),matrix(0,p,1))
      x <- selector*xx + abs(selector-1)*(y-alpha*dfy)
    }
  }
  retOb <- structure(
    list(call = match.call(),
         x = x,
         k = k),
    class = "APG_EN2")
  return(retOb)
}
