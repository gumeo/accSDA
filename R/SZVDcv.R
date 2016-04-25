#' Cross-validation of sparse zero variance discriminant analysis
#'
#' Applies alternating direction methods of multipliers to solve sparse
#' zero variance descriminant analysis.
#'
#' @param X n by p data matrix.
#' @param Y n by K indicator matrix.
#' @param folds number of folds to use in K-fold cross-validation.
#' @param gams Vector of regularization parameters for l1 penalty.
#'        All elements must be greater than zero.
#' @param beta Augmented Lagrangian parameter. Must be greater than zero.
#' @param D Penalty dictionary basis matrix.
#' @param q Desired number of discriminant vectors.
#' @param maxits Number of iterations to runn ADMM algorithm.
#' @param tol Stopping tolerances for ADMM, must have tol$rel and tol$abs.
#' @param ztol Rounding tolerance for truncating entries to 0.
#' @param feat Maximum fraction of nonzero features desired in validation scheme.
#' @param quiet toggles between displaying intermediate statistics.
#' @return \code{SZVDcv} returns an object of \code{\link{class}} "\code{SZVDcv}" 
#'        including a list with the named components \code{DVs} and \code{gambest}.
#'        Where \code{DVs} are the discriminant vectors for the best l1 regularization
#'        parameter and \code{gambest} is the best regularization parameter found
#'        in the cross-validation.
#' @seealso Used by: \code{\link{SZVDcv}}.
#' @details
#' This function will currently solve as a standalone function in accSDA for time comparison.
#' A wrapper function like ASDA will be created to use the functionality of plots and such.
#' Maybe call it ASZDA. For that purpose the individual ZVD function will need to be implemented.
#' @rdname SZVDcv
#' @export SZVDcv
SZVDcv <- function(X, ...) UseMethod("SZVDcv")

#' @return \code{NULL}
#'
#' @rdname SZVDcv
#' @method SZVDcv default
SZVDcv.default <- function(X, Y, folds, gams,  beta,D, q, maxits, tol, ztol, feat, quiet){
  ## Initialization
  
  # Gams is a vector, change into a matrix
  if(is.null(dim(gams))){
    tmpG <- matrix(0,length(gams),q)
    for(i in 1:q){
      tmpG[,i] <- gams
    }
    gams <- tmpG
  }
  
  # Get dimensions of input matrices
  n <- dim(X)[1]
  p <- dim(X)[2]
  K <- dim(Y)[2]
  
  # If n is not divisible by K, duplicate some records for the sake of
  # cross validation.
  pad <- 0
  if(n %% folds > 0){
    pad <- ceiling(n/folds)*folds - n
    
    # Add the duplicates, such that number of data points is
    # divisible by the number of folds
    X <- rbind(X,X[1:pad,])
    Y <- rbind(Y,Y[1:pad,])
  }
  
  # Get the new number of rows
  n <- dim(X)[1]
  
  # Randomly permute rows of X
  prm <- sample(1:n,n,replace=FALSE)
  X <- X[prm,]
  Y <- Y[prm,]
  
  # Make Atrain
  Atrain <- cbind(Y%*%matrix(1:K,K,1),X)
  
  ###
  # Initialization of cross-validation indices
  ###
  
  # Number of validation samples
  nv <- n/folds
  
  # Initial validation indices
  vinds <- 1:nv
  
  # Initial training indices
  tinds <- (nv+1):n
  
  # Number of params to test
  ngam <- dim(gams)[1]
  
  # Validation scores
  scores <- q*p*matrix(1,nrow = folds, ncol = ngam)
  
  # Misclassification rate for each classifier
  mc <- matrix(0,nrow = folds, ncol = ngam)
  
  for(f in 1:folds){
    ## Initialization
    
    # Extract X and Y data
    Xt <- X[tinds,]
    At <- Atrain[tinds,]
    
    # Extract validation data
    Av <- Atrain[vinds,]
    
    # Get dimensions of training matrices.
    nt <- dim(Xt)[1]
    p <- dim(Xt)[2]
    
    # Call ZVD function to solve the unpenalized problem
    get_DVs <- TRUE
    w0 <- ZVD(At, 0, get_DVs)
    
    # Normalize B (divide by the spectral norm)
    if(dim(w0$B)[2] == 1){ # B stored as vector
      w0$B <- w0$B/norm(w0$B, type = "2")
    } else{ # B stored as p by p symmetric matrix
      w0$B <- w0$B + t(w0$B)
      w0$B <- w0$B/norm(w0$B, type = "2")
    }
    
    # Initialize objective matrix
    if(dim(w0$B)[2] == 1){ # B stored as vector
      B0 <- w0$B
    } else{ # B stored as p by p symmetric matrix
      B0 <- crossprod(w0$N,w0$B)%*%w0$N
      B0 <- (B0 + t(B0))/2
    }
    
    # Initialize the nullspace matrix
    N0 <- w0$N
    
    # Number of gammas
    num_gammas <- dim(gams)[1]
    
    # Validation loop
    if(!quiet){
      print("-------------------------------------------")
      print(paste("Fold number:",f))
      print("-------------------------------------------")
    }
    
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Loop through potential regularization parameters.
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    for(ll in 1:num_gammas){
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      # Initialization.
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      # Initialize B and N
      B <- B0
      N <- N0
      
      # Initialize DVs
      DVs <- matrix(0,p,q)
      
      # Set x0 to be the unpenalized zero-variance discriminant vector in Null(W0)
      if(dim(B0)[2] == 1){ #B0 vector
        # Compute t(DN)%*%(mu1-mu2)
        w <- crossprod(N0,crossprod(D,B0))
        
        # Use normalized w as initial x
        x0 <- w/norm(w,type = "2")
      } else{ #B0 matrix
        x0 <- crossprod(N0,crossprod(D,w0$dvs[,1,drop=FALSE]))
      }
      
      # y is the unpenalized solution in the original space.
      # z is the all-zeros vector in the original space.
      
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      # Call Alternating Direction Method to solve SZVD problem.
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      for(j in 1:q){
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Initialization.
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        sols0 <- list(x = x0,
                      y = w0$dvs[,j],
                      z = matrix(0,p,1))
        quietADMM <- TRUE
        
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Call ADMM solver.
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++
        s <- matrix(1,p,1)
        retObj <- SZVD_ADMM(B, N, D, sols0, s, gams[ll,j], beta, tol, maxits, quietADMM)
        tmpx <- retObj$x
        
        # Extract j-th discriminant vector
        DVs[,j] <- D%*%N%*%tmpx
        
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Update N and B for the newly found DV.
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        if(j < q){
          # Project columns of N onto orthogonal complement of Nx
          
          x <- DVs[,j]
          x <- x/norm(x, type = "2")
          
          # Project N into orthogonal complement of span(x)
          Ntmp <- N - x%*%(crossprod(x,N))
          
          # Call QR factorization to extract orthonormal basis for span (Ntmp)
          qrTmp <- qr(Ntmp)
          Q <- qr.Q(qrTmp)
          R <- qr.R(qrTmp)
          
          # Extract nonzero rows of R to get columns of Q to use as new N.
          R_rows <- (abs(diag(R)) > 1e-6)
          
          # Use nontrivial columns of Q as updated N.
          N <- Q[, R_rows]
          
          # Update B0 according to the new basis N.
          B <- crossprod(N, w0$B) %*% N
          B <- 0.5*(B+t(B))


          # Update initial solutions in x direction by projecting next unpenalized ZVD vector.
          x0 <- crossprod(N, crossprod(D, w0$dvs[,j+1]))
        } 
        
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Get performance scores on the validation set.
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Call test_ZVD to get predictions, etc.
        
        statsObj <- test_ZVD(DVs, Av, w0$means, w0$mu, FALSE)
        stats <- statsObj$stats
        
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Validation scores.
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        # Round small entries to zero
        DVs <- DVs * (ceiling(abs(DVs)-ztol))
        
        # if fraction nonzero features less than feat.
        if( 1 <= sum(DVs != 0) & sum(DVs != 0) <= q*p*feat){
          # Use misclassification rate as validation score.
          scores[f,ll] <- stats$mc
        } else if(sum(DVs != 0) < 0.5){ 
          # Trivial solution
          scores[f,ll] <- q*p # Disqualify with largest possible value
        } else{
          # Solution is not sparse enough, use most sparse as measure of quality instead.
          scores[f,ll] <- sum(DVs != 0)
        }
        
        # Display iteration stats.
        if(!quiet){
          print(paste("f:", f, "| ll:", ll, "| lam:", lams[ll], "| feat:", 
                      sum(DVs != 0)/(q*p), "| mc:", stats$mc, "| score:", scores[f,ll]))
        }
      }
    } # End ll over validation parameters
    #--------------------------------------------
    # Update training/validation split
    #--------------------------------------------
    # Extract new validation indices
    tmp <- tinds[1:nv]
    
    if(nv+1 > nt){
      # Special case for 2-fold CV
      tinds <- vinds
      vinds <- tmp
    } else{
      tinds <- c(tinds[(nv+1):nt],vinds)
      
      # Update validation indices
      vinds <- tmp
    }
  } # End folds
  
  ###
  # Find the best solution
  ###
  
  # Average CV scores
  avg_score <- colMeans(scores)
  
  # Choose lambda with best average score
  gbest <- which.min(avg_score)
  
  gambest <- gams[gbest]
  
  ###
  # Solve with lambda = lam(lbest)
  ###
  print(paste("Finished Training: lam =", gambest))
  
  # Use the full training set to obtain parameters
  Xt <- X[1:(n-pad),]
  Yt <- Y[1:(n-pad),]
  Atrain <- cbind(Yt%*%matrix(1:K,K,1),Xt)
  # Loop until nontrivial solution is found
  trivsol <- TRUE
  while(trivsol){
    szvdObj <- SZVD(Atrain, gambest, D, FALSE, FALSE, tol, maxits, beta, quietADMM)
    DVs <- szvdObj$DVs
    
    # Round small entried to zero
    DVs <- DVs*(ceiling(abs(DVs)-ztol))
    
    # Check for trivial solution
    if(sum(DVs != 0) == 0){
      # If trivial solution, update gbest by one and update gambest.
      gbest <- gbest + 1
      gambest <- gams[gbest,]
    } else{
      trivsol <- FALSE
    }
  }
  
  # Create an object of class SDAPcv to return, might add more to it later
  retOb <- structure(
    list(call = match.call(),
         DVs = DVs,
         gambest = gambest),
    class = "SZVDcv")
  
  return(retOb)
}

#' @export
SZVDcv.matrix <- function(X, ...){
  res <- SZVDcv.default(X, ...)
  cl <- match.call()
  cl[[1L]] <- as.name("SZVDcv")
  res$call <- cl
  res
}