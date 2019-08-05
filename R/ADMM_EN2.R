#' ADMM on l1 regularized quadratic program
#'
#' Applies Alternating Direction Method of Multipliers to the l1-regularized quadratic program
#' \deqn{f(\mathbf{x}) + g(\mathbf{x}) = \frac{1}{2}\mathbf{x}^TA\mathbf{x} - d^T\mathbf{x} + \lambda |\mathbf{x}|_1}{f(x) + g(x) = 0.5*x^T*A*x - d^T*x + lambda*|x|_l1}
#'
#' @param R Upper triangular matrix in Chol decomp \eqn{\mu I + A = R^T R}{mu*I + A = R'*R}.
#' @param d nx1 dimensional column vector.
#' @param lam Regularization parameter for l1 penalty, must be greater than zero.
#' @param mu Augmented Lagrangian penalty parameter, must be greater than zero.
#' @param alpha Step length.
#' @param maxits Number of iterations to run
#' @param tol Vector of stopping tolerances, first value is absolute, second is relative tolerance.
#' @param quiet Logical controlling display of intermediate statistics.
#' @param selector Vector to choose which parameters in the discriminant vector will be used to calculate the
#'                 regularization terms. The size of the vector must be *p* the number of predictors. The
#'                 default value is a vector of all ones. This is currently only used for ordinal classification.
#' @return \code{ADMM_EN2} returns an object of \code{\link{class}} "\code{ADMM_EN2}" including a list
#' with the following named components
#'
#' \describe{
#'   \item{\code{call}}{The matched call.}
#'   \item{\code{x}}{Found solution.}
#'   \item{\code{y}}{Dual solution.}
#'   \item{\code{z}}{Slack variables.}
#'   \item{\code{k}}{Number of iterations used.}
#' }
#' @seealso Used by: \code{\link{SDAD}} and the \code{SDADcv} cross-validation version.
#' @details
#' This function is used by other functions and should only be called explicitly for
#' debugging purposes.
#' @keywords internal
ADMM_EN2 <- function(R, d, x0, lam, mu, maxits, tol, quiet, selector = rep(1,dim(x)[1])){
  ###
  # Initialization
  ###
  x <- x0
  y <- x0
  p <- dim(x)[1]
  z <- matrix(0,p,1)

  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Outer loop: Repeat until converged or max # of iterations reached.
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  for(k in 0:maxits){
    ###
    # Update x
    ###
    # (mu I + A)x = d + mu*y - z
    b <- d + mu*y - z
    #Rx <- solve(t(R),b)
    #x <- solve(R,Rx)
    Rx <- forwardsolve(t(R),b)
    x <- backsolve(R,Rx)

    ###
    # Update y
    ###
    # Update using soft-thresholding
    yold <- y
    tmp <- x + z/mu
    yy <- sign(tmp)*pmax(abs(tmp)-(lam/mu)*matrix(1,p,1),matrix(0,p,1))
    y <-  selector*yy + abs(selector-1)*(tmp)

    ###
    # Update z
    ###
    z <- z + mu*(x-y)

    ###
    # Check convergence
    ###
    # Primal constraint violation

    # Primal residual
    r <- x - y

    # l2 norm of residual
    dr <- norm(r, type = "2")

    # Dual constraint violation

    # Dual residual
    s <- mu*(y - yold)

    # l2 norm of the residual
    ds <- norm(s, type = "2")

    ###
    # Check if stopping criteria is satisfied
    ###

    # Compute absolute and relative tolerances
    ep <- sqrt(p)*tol[1] + tol[2]*max(norm(x, type = "2"), norm(y, type = "2"))
    es <- sqrt(p)*tol[1] + tol[2]*norm(y, type = "2")

    # Display iteration stats
    if(!quiet){
      print(paste("it = ", k, ", primal_viol = ",
                  dr-ep, ", dual_viol = ", ds-es,
                  ", norm_y = ",
                  max(norm(x, type = "2"), norm(y, type = "2")), sep=""))
    }

    # Check if the residual norms are less than the given tolerance
    if(dr < ep & ds < es){
      break # Convergence
    }

  }
  retOb <- structure(
    list(call = match.call(),
         x = x,
         y = y,
         z = z,
         k = k),
    class = "ADMM_EN2")
  return(retOb)
}
