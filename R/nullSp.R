#' Finding null space of linear operator
#'
#' Finds the null space of a linear operator A in \eqn{R^{n \times m}}{R^(n by m)}.
#' The null space is given as a matrix, where the columns form an orthonormal basis
#' for the nullspace. This function emulates the null function in matlab, it
#' works exactly the same, but the basis vectors may be different, i.e. rotated.
#'
#' @param A m by n matrix
#' @return \code{nullSp} returns a matrix whose columns span the nullspace of A.
#' @seealso Alternative \code{\link[MASS]{Null}} function in MASS package.
#' @details
#' This function is used by other functions and should only be called explicitly for
#' debugging purposes.
#' @keywords internal
nullSp <- function(A){
  m <- dim(A)[1]
  n <- dim(A)[2]
  val <- svd(A, nv = n)
  V <- val$v
  S <- val$d
  if(m > 1){
    s <- S
  } else if(m == 1){
    s <- S[1]
  } else{
    s <- 0
  }
  # This is the matlab code for the line below
  # here we hardcode that we use double precision
  #tol = max(m,n) * max(s) * eps(class(A));
  tol <- max(m,n)*max(s)*.Machine$double.eps
  r <- sum(s > tol)
  # Return empty matrix if nothing in nullspace
  if(r+1>n){
    N <- matrix(0,dim(V)[1],0)
  } else{
    N <- V[,(r+1):n, drop = FALSE] # The null space
  }
  return(N)
}
