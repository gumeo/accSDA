#' Normalize training data
#'
#' Normalize test data using output from the \code{normalize()} of the training data
#'
#' @param Xtst a matrix with the test data with observations down the rows and variables in the columns.
#' @param Xn List with the output from normalize(Xtr) of the training data.
#' @return \code{normalizetest} returns the normalized test data \code{Xtst}
#' @author Line Clemmensen
#' @references Clemmensen, L., Hastie, T. and Ersboell, K. (2008)
#' "Sparse discriminant analysis", Technical report, IMM, Technical University of Denmark
#' @details
#' This function can e.g. be used for the test data in the \code{\link{predict.ASDA}} function.
#' @seealso \code{\link{normalize}}, \code{\link{predict.ASDA}}, \code{\link{ASDA}}
#' @examples
#' ## Data
#' Xtr<-matrix(sample(seq(3),12,replace=TRUE),nrow=3)
#' Xtst<-matrix(sample(seq(3),12,replace=TRUE),nrow=3)
#'
#' ## Normalize training data
#' Nm<-normalize(Xtr)
#'
#' ## Normalize test data
#' Xtst<-normalizetest(Xtst,Nm)
#'
#' @rdname normalizetest
#' @export normalizetest
normalizetest <- function(Xtst,Xn){
  ntst <- dim(Xtst)[1]
  p <- dim(Xtst)[2]
  m <- rep(Xn$mx,ntst)
  dim(m) <- c(p,ntst)
  Xtst <- Xtst-t(m)
  v <- rep(Xn$vx,ntst)
  dim(v) <- c(p,ntst)
  v <- t(v)
  Xtst <- Xtst[,Xn$Id]/v[1:ntst,Xn$Id]
}

