#' ADMM on l1 regularized quadratic program
#'
#' Applies Alternating Direction Method of Multipliers to the l1-regularized quadratic program
#' \deqn{f(\mathbf{x}) + g(\mathbf{x}) = \frac{1}{2}\mathbf{x}^TA\mathbf{x} - d^T\mathbf{x} + \lambda |\mathbf{x}|_1}{f(x) + g(x) = 0.5*x^T*A*x - d^T*x + lambda*|x|_l1}
#'
#' @param Ainv Diagonal of \eqn{A^{-1}}{A^{-1}} term in SMW formula, where A is an n by n
#' positive definite coefficient matrix.
#' @param V Matrix from SMW formula.
#' @param R Upper triangular matrix in Cholesky decomposition of \eqn{I + UA^{-1}V}{I + U*Ainv*V}.
#' @param d nx1 dimensional column vector.
#' @param lam Regularization parameter for l1 penalty, must be greater than zero.
#' @param mu Augmented Lagrangian penalty parameter, must be greater than zero.
#' @param alpha Step length.
#' @param maxits Number of iterations to run
#' @param tol Vector of stopping tolerances, first value is absolute, second is relative tolerance.
#' @param quiet Logical controlling display of intermediate statistics.
#' @return \code{ADMM_EN_SMW} returns an object of \code{\link{class}} "\code{ADMM_EN_SMW}" including a list
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
ADMM_EN_SMW <- function(Ainv, V,R, d, x0, lam, mu, maxits, tol, quiet){
  ###
  # Initialization
  ###

  x <- x0
  y <- x0
  p <- dim(x)[1]
  z <- matrix(0,p,1)
  n <- dim(V)[1]

  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Outer loop: Repeat until converged or max # of iterations reached.
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  for(k in 0:maxits){
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Update x using SMW applied to
    # (mu I + gam*Om + X'X)x = d + mu*y - z.
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    # RHS coefficient vectors.
    b <- d + mu*y - z
    btmp <- V%*%(Ainv*b)/n # Ainv is a vector representing a diagonal matrix

    # Apply SMW to get x
    #x <- Ainv*b - 2*Ainv*(t(V)%*%(solve(R,solve(t(R),btmp))))
    x <- Ainv*b - 2*Ainv*(t(V)%*%(backsolve(R,forwardsolve(t(R),btmp))))

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Update y.
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Update y using soft-thresholding.
    yold <- y
    tmp <- x + z/mu
    y <- sign(tmp)*pmax(abs(tmp) - lam*matrix(1,p,1),matrix(0,p,1))

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Update z.
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    z  <- z + mu*(x-y)

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Check for convergence.
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ### Primal constraint violation.

    # Primal residual.
    r <- x - y

    # l2 norm of the residual.
    dr <- norm(r, type = "2")

    ### Dual constraint violation.

    # Dual residual.
    s <- mu*(y - yold)

    # l2 norm of the residual.
    ds <- norm(s, type = "2")

    ###  Check if the stopping criteria are satisfied.

    # Compute absolute and relative tolerances.
    ep = sqrt(p)*tol[1] + tol[2]*max(norm(x, type = "2"), norm(y, type = "2"))
    es = sqrt(p)*tol[1] + tol[2]*norm(y, type = "2")

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
    class = "ADMM_EN_SMW")
  return(retOb)
}
