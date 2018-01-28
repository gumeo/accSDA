#' Accelerated Proximal Gradient on l1 regularized quadratic program with backtracking
#'
#' Applies accelerated proximal gradient (with backtracking) algorithm to the l1-regularized quadratic program
#' \deqn{f(\mathbf{x}) + g(\mathbf{x}) = \frac{1}{2}\mathbf{x}^TA\mathbf{x} - d^T\mathbf{x} + \lambda |\mathbf{x}|_1}{f(x) + g(x) = 0.5*x^T*A*x - d^T*x + lambda*|x|_l1}
#'
#' @param A p by p positive definite coefficient matrix
#' \deqn{A = (\gamma Om + X^T X/n)}{A = (gamma Om + X^T X/n)}.
#' @param Xt Same as X above, we need it to make calculations faster.
#' @param Om Same reason as for the above parameter.
#' @param gamma l2 regularizing parameter.
#' @param d nx1 dimensional column vector.
#' @param lam Regularization parameter for l1 penalty, must be greater than zero.
#' @param L Initial value of backtracking Lipshitz constant.
#' @param eta Backtracking scaling parameter.
#' @param maxits Number of iterations to run
#' @param tol Stopping tolerance for proximal gradient algorithm.
#' @return \code{prox_ENbt} returns an object of \code{\link{class}} "\code{prox_ENbt}" including a list
#' with the following named components
#'
#' \describe{
#'   \item{\code{call}}{The matched call.}
#'   \item{\code{x}}{Found solution.}
#'   \item{\code{k}}{Number of iterations used.}
#' }
#' @seealso Used by: \code{\link{SDAP}} and the \code{SDAPcv} cross-validation version.
#' @details
#' This function is used by other functions and should only be called explicitly for
#' debugging purposes.
#' @keywords internal
prox_ENbt <- function(A, Xt, Om, gamma,  d, x0, lam, L, eta, maxits, tol){
  # Make sure these are not a matrix with
  # one element
  origL <- L
  lam <- as.numeric(lam)
  #alpha <- as.numeric(alpha)
  ifDiag <- FALSE
  if(norm(diag(diag(Om))-Om, type = "F") < 1e-15){
    ifDiag <- TRUE
  }else{
    # Factorize Omega
    R <- chol(gamma*Om,pivot=TRUE)
  }

  ###
  # Initialization
  ###
  # Initial solution
  x <- x0

  # Get number of components of x,d, rows/cols of A
  n <- dim(x)[1]

  oneMat <- matrix(1,n,1)
  zeroMat <- matrix(0,n,1)

  ###
  # Outer loop: Repeat until converged or max # of iterations reached.
  ###
  for(k in 0:maxits){
    # Compute gradient of differentiable part (f(x) = 0.5*x'*A*x - d'*x)
    df <- A%*%x - d

    ###
    # Compute disagreement between df and lam*sign(x) on supp(x).
    ###
    # Initialize error vector
    err <- matrix(0,nrow = n, ncol = 1)
    # Initialize cardinality of support
    card <- 0

    # For each i, update error if i in the support
    for(i in 1:n){
      if(abs(x[i]) > 1e-12){
        # Update cardinality
        card <- card + 1

        # Update error vector
        err[i] <- -df[i]-lam*sign(x[i])
      }
    }

    ###
    # Check stopping criteria -df(x) in subdiff g(x).
    # Need inf(df) < lam + tol, inf(err) < tol.
    ###
    if(max(norm(df,type="I")-lam, norm(err,type="I")) < tol*n){
      # CONVERGED!
      return(structure(
        list(call = match.call(),
             x = x,
             k = k,
             L = L),
        class = "prox_ENbt"))
    } else{
      #--------------------------------------------------------------------------------------
      # Backtracking line search update
      #--------------------------------------------------------------------------------------
      #L <- origL
      alpha <- 1/L # step length
      # Evaluate proximal gradient
      pL <- sign(x-alpha*df)*pmax(abs(x-alpha*df) - lam*alpha*oneMat,zeroMat)
      pTilde <- (pL-x)
      if(ifDiag == TRUE){
        gap <- as.numeric((1/2)*(sum(pTilde*pTilde*(L*rep(1,n)-gamma*diag(Om)))-norm(Xt%*%pTilde,type='2')^2))
      }else{
        gap <- as.numeric((1/2)*(sum(pTilde*pTilde*L)-norm(R%*%pTilde,type='2')^2-norm(Xt%*%pTilde,type='2')^2))
        #gap <- as.numeric((1/2)*t(pTilde)%*%(L*diag(n)-A)%*%(pTilde))
      }

      # backtrack
      while(gap < -tol){
        L <- eta*L
        alpha <- 1/L # step length
        # Evaluate proximal gradient
        pL <- sign(x-alpha*df)*pmax(abs(x-alpha*df) - lam*alpha*oneMat,zeroMat)
        pTilde <- (pL-x)
        if(ifDiag == TRUE){
          gap <- as.numeric((1/2)*(sum(pTilde*pTilde*(L*rep(1,n)-gamma*diag(Om)))-norm(Xt%*%pTilde,type='2')^2))
        }else{
          gap <- as.numeric((1/2)*(sum(pTilde*pTilde*L)-norm(R%*%pTilde,type='2')^2-norm(Xt%*%pTilde,type='2')^2))
          #gap <- as.numeric((1/2)*t(pTilde)%*%(L*diag(n)-A)%*%(pTilde))
        }
      }
      x <- pL
    }
  }
  retOb <- structure(
    list(call = match.call(),
         x = x,
         k = k,
         L = L),
    class = "prox_ENbt")
  return(retOb)
}
