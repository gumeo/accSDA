#' Accelerated Proximal Gradient (with backtracking) on l1 regularized quadratic program
#'
#' Applies accelerated proximal gradient algorithm (with backtracking) to the l1-regularized quadratic program
#' \deqn{f(\mathbf{x}) + g(\mathbf{x}) = \frac{1}{2}\mathbf{x}^TA\mathbf{x} - d^T\mathbf{x} + \lambda |\mathbf{x}|_1}{f(x) + g(x) = 0.5*x^T*A*x - d^T*x + lambda*|x|_l1}
#'
#' @param A p by p positive definite coefficient matrix
#' \deqn{A = (\gamma Om + X^T X/n)}{A = (gamma Om + X^T X/n)}.
#' @param Xt Same as X above, we need it to make calculations faster.
#' @param Om Same reason as for the above parameter.
#' @param gamma l2 regularizing parameter.
#' @param d nx1 dimensional column vector.
#' @param lam Regularization parameter for l1 penalty, must be greater than zero.
#' @param L Initial vlaue of the backtracking Lipshitz constant.
#' @param eta Backtracking scaling parameter.
#' @param maxits Number of iterations to run
#' @param tol Stopping tolerance for proximal gradient algorithm.
#' @param selector Vector to choose which parameters in the discriminant vector will be used to calculate the
#'                 regularization terms. The size of the vector must be *p* the number of predictors. The
#'                 default value is a vector of all ones. This is currently only used for ordinal classification.
#' @return \code{APG_EN2bt} returns an object of \code{\link{class}} "\code{APG_EN2bt}" including a list
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
APG_EN2bt <- function(A, Xt, Om, gamma, d, x0, lam, L, eta,  maxits, tol, selector = rep(1,dim(x0)[1])){
  ifDiag <- FALSE
  if(norm(diag(diag(Om))-Om, type = "F") < 1e-15){
    ifDiag <- TRUE
  }else{
    # Factorize Omega
    R <- chol(gamma*Om,pivot=TRUE)
  }
  origL <- L
  ###
  # Initialization
  ###

  x <- x0
  y <- x

  # Get number of components of x
  p <- dim(x0)[1]

  # initial momentum coefficient
  t <- 1
  #told <- 1

  # Objective function and gradient
  if(A$flag == 1){
    #f <- function(x){
    #  as.numeric(t(x)%*%(matrix(A$gom,length(A$gom),1)*x) + (1/A$n)*norm(A$X%*%x)^2 + t(d)%*%x)
    #}
    df <- function(x){
      #2*(matrix(A$gom,length(A$gom),1)*x + crossprod(A$X,A$X%*%(x/A$n))) - d
      2*(matrix(A$gom,length(A$gom),1)*x + crossprod(A$X,A$X%*%(x/A$n))) - d
    }
  }else{
    #f <- function(x){
    #  as.numeric(0.5*t(x)%*%A$A%*%x + t(d)%*%x)
    #}
    df <- function(x){
      A$A%*%x - d
    }
  }
  oneMat <- matrix(1,p,1)
  zeroMat <- matrix(0,p,1)
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
      # Update gradient via backtracking
      ###
      #L <- origL
      alpha <- 1/L # step length
      dfy <- df(y)
      pLyy <- sign(y-alpha*dfy)*pmax(abs(y-alpha*dfy) - lam*alpha*oneMat,zeroMat)
      pLy <- selector*pLyy + abs(selector-1)*(y-alpha*dfy)

      pTilde <- (pLy-y)
      if(ifDiag == TRUE){
        #QminusF = 1/2*L*norm(dy,2)^2 - 1/A.n*norm(A.X*dy)^2 - dy'*((A.gom).*dy);
        QminusF <- as.numeric((0.5)*L*norm(pTilde, type = "2")^2 - (1/A$n)*norm(A$X%*%pTilde, type="2")^2 - sum((pTilde)*(A$gom*pTilde)))
      }else{
        #QminusF = 1/2*(L*norm(dy, 2)^2 - dy'*A.A*dy);
        QminusF <- as.numeric((1/2)*(L*norm(pTilde, type="2")^2 - t(pTilde)%*%A$A%*%pTilde))
      }

      #QminusF <- as.numeric((1/2)*t(pLy-y)%*%(L*diag(p)-A$A)%*%(pLy-y))
      while(QminusF < -tol){
        L <- eta*L
        alpha <- 1/L # step length
        pLyy <- sign(y-alpha*dfy)*pmax(abs(y-alpha*dfy) - lam*alpha*oneMat,zeroMat)
        pLy <- selector*pLyy + abs(selector-1)*(y-alpha*dfy)
        pTilde <- (pLy-y)
        if(ifDiag == TRUE){
          QminusF <- as.numeric((0.5)*L*norm(pTilde, type = "2")^2 - (1/A$n)*norm(A$Xt%*%pTilde, type="2")^2 - sum((pTilde)*(A$gom*pTilde)))
        }else{
          QminusF <- as.numeric((1/2)*(L*norm(pTilde, type="2")^2 - t(pTilde)%*%A$A%*%pTilde))
        }

        #QminusF <- as.numeric((1/2)*t(pLy-y)%*%(L*diag(p)-A$A)%*%(pLy-y))
      }

      # Update based on backtracked APG solution
      xold <- x
      x <- pLy

      # Extrapolation factor
      told <- t
      t <- (1 + sqrt(1+4*told^2))/2

      # Extrapolate using last two iterates
      y <- x+(told-1)/t*(x-xold)
    }
  }
  retOb <- structure(
    list(call = match.call(),
         x = x,
         k = k,
         L = L),
    class = "APG_EN2bt")
  return(retOb)
}
