#' Normalize training data
#'
#' Normalize a vector or matrix to zero mean and unit length columns.
#'
#' @param X a matrix with the training data with observations down the rows and variables in the columns.
#' @return \code{normalize} Returns a list with the following attributes:
#'     \describe{
#'       \item{Xc}{The normalized data}
#'       \item{mx}{Mean of columns of \code{X}.}
#'       \item{vx}{Length of columns of \code{X}.}
#'       \item{Id}{Logical vector indicating which variables are
#'       included in X. If some of the columns have zero length they are omitted}
#'     }
#' @author Line Clemmensen
#' @references Clemmensen, L., Hastie, T. and Ersboell, K. (2008)
#' "Sparse discriminant analysis", Technical report, IMM, Technical University of Denmark
#' @details
#' This function can e.g. be used for the training data in the \code{ASDA} function.
#' @seealso \code{\link{normalizetest}}, \code{\link{predict.ASDA}}, \code{\link{ASDA}}
#' @examples
#' ## Data
#' X<-matrix(sample(seq(3),12,replace=TRUE),nrow=3)
#'
#' ## Normalize data
#' Nm<-normalize(X)
#' print(Nm$Xc)
#'
#' ## See if any variables have been removed
#' which(!Nm$Id)
#' @rdname normalize
#' @export normalize
normalize <- function(X){
  n <- dim(X)[1]
  p <- dim(X)[2]
  mx <- apply(X,2,mean)
  m <- rep(mx,n)
  dim(m) <- c(p,n)
  X <- X-t(m)
  vx <- sqrt(apply(X^2,2,sum))
  Id <- vx!=0
  v <- rep(vx,n)
  dim(v) <- c(p,n)
  v <- t(v)
  X <- X[,Id]/v[,Id]
  list(Xc=X,mx=mx,vx=vx,Id=Id)
}


