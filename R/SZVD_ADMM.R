#' Alternating Direction Method of Multipliers for SZVD
#'
#' Iteratively solves the problem
#' \deqn{\text{min}(-1/2*x^TB^Tx + \gamma p(y): ||x||_2 \leq 1, DNx = y)}{min(-1/2*x^TB^Tx + gamma p(y): l2(x) <= 1, DNx = y)}
#'
#' @param B Between class covariance matrix for objective (in space defined by N).
#' @param N basis matrix for null space of covariance matrix W.
#' @param D penalty dictionary/basis.
#' @param sols0 initial solutions sols0$x, sols0$y, sols0$z
#' @param pen_scal penalty scaling term.
#' @param gamma l1 regularization parameter
#' @param beta penalty term controlling the splitting constraint.
#' @param tol tol$abs = absolute error, tol$rel = relative error to be
#'        achieved to declare convergence of the algorithm.
#' @param maxits maximum number of iterations of the algorithm to run.
#' @param quiet toggles between displaying intermediate statistics.
#' @return \code{SZVD_ADMM} returns an object of \code{\link{class}} "\code{SZVD_ADMM}" including a list
#' with the following named components
#'
#' \describe{
#'   \item{\code{x,y,z}}{Iterates at termination.}
#'   \item{\code{its}}{Number of iterations required to converge.}
#'   \item{\code{errtol}}{Stopping error bound at termination}
#' }
#' @seealso Used by: \code{\link{SZVDcv}}.
#' @details
#' This function is used by other functions and should only be called explicitly for
#' debugging purposes.
#' @keywords internal
SZVD_ADMM <- function(B,  N, D, sols0, pen_scal, gamma, beta, tol, maxits, quiet = TRUE){
  #====================================================================
  # Precomputes quantities that will be used repeatedly by the algorithm.
  #====================================================================
  # Dimension of decision variable
  p <- dim(D)[1]

  # Compute D%*%N
  if(all(D == diag(p))){
    DN <- N
  } else{
    DN <- D%*%N
  }

  #====================================================================
  # Initialize x solution and constants for x update.
  #====================================================================
  if(dim(B)[2] == 1){ # Case where K = 2
    # Compute (DN)'*(mu1-mu2)
    w <- crossprod(DN,B)

    # Compute initial x if omitted.
    if (missing(sols0)){
      x = w/norm(w,'f')
    }
    else{
      x = sols0$x
    }

    # Constant for x update step
    Xup <- beta - crossprod(w)
    Xup <- as.numeric(Xup[1,1])
  } else{ # K > 2 case
    # Dimension of the null-space of W.
    l <- dim(N)[2]

    # Compute initial x if omitted.
    if (missing(sols0)){
      x = eigen((B+t(B))/2)$vectors[,1]
    }
    else{
      x = sols0$x
    }

    # Take Cholesky of beta I - B (for use in update of x)
    L = t(chol((beta*diag(l) - B)));
  }

  #====================================================================
  # Initialize decision variables y and z.
  #====================================================================

  # If omitted, set initial y equal to unpenalized ZVD, and z = 0 to start.
  if (missing(sols0)){
    y = DN %*%x
    z = as.matrix(rep(0,p))
  }
  else{
    y = sols0$y
    z = sols0$z
  }

  #====================================================================
  # Call the algorithm.
  #====================================================================

  for(iter in 1:maxits){
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Step 1: Perform shrinkage to update y_k+1.
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    # Save previous iterate
    yold <- y

    # Call soft-thresholding
    y <- vec_shrink(beta*DN%*%x + z, gamma*pen_scal)

    # Normalize y (if necessary).
    tmp = max(c(0, norm(as.matrix(y), 'f') - beta ))
    y = y/(beta + tmp);

    # Truncate complex part (if have +0i terms appearing.)
    y = Re(y)

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Step 2: Update x_k+1 by solving
    # x_k+1 = argmin { -x'*A*x + beta/2 l2(x - y_k+1 + z_k)^2}
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    # Compute RHS
    b <- crossprod(DN,beta*y-z)

    if(dim(B)[2] == 1){ # K = 2
      # Update using Sherman-Morrison formula
      x <- (b + ((t(b)%*%w)[1,1]/Xup)*w)/beta
    } else{ # K > 2
      # Update by solving the system L*t(L)*x = b
      btmp <- forwardsolve(L,b)
      x <- backsolve(t(L),btmp)
    }

    # Truncate complex part (if have +0i terms appearing.)
    x = Re(x)

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #  Step 3: Update the Lagrange multipliers
    # (according to the formula z_k+1 = z_k + beta*(N*x_k+1 - y_k+1) ).
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    zold <- z
    z <- Re(z + beta*(DN%*%x - y))

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Check stopping criteria.
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    #----------------------------------------------------------------
    # Primal constraint violation.
    #----------------------------------------------------------------
    # Primal residual.
    r <- DN%*%x - y
    # l2 norm of the residual
    dr <- norm(as.matrix(r), type = "f")

    #----------------------------------------------------------------
    # Dual constraint violation.
    #----------------------------------------------------------------
    # Dual residual.
    s <- beta*(y-yold)
    # l2 norm of the residual
    ds <- norm(as.matrix(s), type = "f")

    #----------------------------------------------------------------
    # Check if the stopping criteria are satisfied.
    #----------------------------------------------------------------

    # Compute absolute and relative tolerances
    ep <- sqrt(p)*tol$abs + tol$rel*max(norm(as.matrix(x), type = "f"), norm(as.matrix(y), type = "f"))
    es <- sqrt(p)*tol$abs + tol$rel*norm(as.matrix(y), type = "f")

    # Display current iteration stats.
    if (quiet==FALSE){
      print(sprintf("it = %g, dr = %e, dr-ep = %e, ds = %e, ds-es = %e, normy = %e",
                    iter, dr, dr-ep, ds, ds-es, norm(as.matrix(y), 'F')), quote=F)
      cat("\n")
    }

    # Check if the residual norms are less than the given tolerance
    if (max(c(dr - ep, ds - es, 0), na.rm=TRUE)==0 & iter > 10){
      break # The algorithm has converged.
    }
  }

  #====================================================================
  # Output results.
  #====================================================================

  return(structure(
    list(call = match.call(),
         x = x,
         y = y,
         z = z,
         its = iter,
         errtol = min(c(ep,es))),
    class = "SZVD_ADMM"))
}

#' Softmax for SZVD ADMM iterations
#'
#' Applies softmax the soft thresholding shrinkage operator to v with tolerance a.
#' That is, output is the vector with entries with absolute value v_i - a if
#' |v_i| > a and zero otherwise, with sign pattern matching that of v.
#'
#' @param v Vector to be thresholded.
#' @param a Vector of tolerances.
#' @return thresholded v vector.
#'
#' \describe{
#'   \item{\code{x,y,z}}{Iterates at termination.}
#'   \item{\code{its}}{Number of iterations required to converge.}
#'   \item{\code{errtol}}{Stopping error bound at termination}
#' }
#' @seealso Used by: \code{\link{SZVD_ADMM}}.
#' @details
#' This function is used by other functions and should only be called explicitly for
#' debugging purposes.
#' @keywords internal
vec_shrink <- function(v,a){
  return(sign(v)*pmax(abs(v)-a, rep(0, length(v))))
}
