#' Cross-validation of sparse zero variance discriminant analysis
#'
#' Applies alternating direction methods of multipliers to solve sparse
#' zero variance descriminant analysis.
#'
#' @param X n by p data matrix.
#' @param Y n by K indicator matrix.
#' @param folds number of folds to use in K-fold cross-validation.
#' @param gams Number of regularly spaced regularization parameters to try in [0,1]*max_gamma.
#'        See details for how max_gamma is computed in the function.
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
#' Add how max_gamma is calculated from the ZVD solution.
#' This function might require a wrapper similar to ASDA.
#' @rdname SZVD_kFold_cv
#' @export SZVD_kFold_cv
SZVD_kFold_cv <- function(X, ...) UseMethod("SZVD_kFold_cv",X)

#' @return \code{NULL}
#'
#' @rdname SZVD_kFold_cv
#' @method SZVD_kFold_cv default
SZVD_kFold_cv.default <- function(X, Y, folds, gams,  beta,D, q, maxits, tol, ztol, feat, quiet){

  # Try to put everything in the style of the old function to make the code for experiments more convenient.
  # Also, let's choose the best parameter in the same way as for the other methods, meaning that we need
  # a feats parameter, setting a threshold for the maximum fraction of nonzero values in the discriminant vectors.
  ## Initialization

  # Get dimensions of input matrices
  n <- dim(X)[1]
  p <- dim(X)[2]
  K <- dim(Y)[2]
  num_gammas <- gams

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
  scores <- q*p*matrix(1,nrow = folds, ncol = num_gammas)

  # Misclassification rate for each classifier
  mc <- matrix(0,nrow = folds, ncol = num_gammas)

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
    w0 <- ZVD(At, scaling = FALSE, get_DVs = TRUE)

    # Extract scaling vector for weighted l1 penalty and diagonal penalty matrix.
    s = sqrt(diag(w0$W))
    w0$s = s

    # If dictionary D missing, use the identity matri.x
    if (missing(D)){
      D = diag(p)
    }

    ## Normalize B (divide by the spectral norm)
    if (dim(w0$B)[2]==1){
      w0$B = w0$B/norm(w0$B, type='f')
    }  else{
      w0$B = (w0$B + t(w0$B))/eigen((w0$B + t(w0$B)), symmetric=TRUE, only.values=TRUE)$values[1]
    }

    # Compute ratio of max gen eigenvalue and l1 norm of the first ZVD to get "bound" on gamma.
    if (dim(w0$B)[2]==1){
      max_gamma =  (t(w0$dvs)%*%w0$B)^2/sum(abs(s*(D %*% w0$dvs)))
    }else{
      max_gamma = apply(w0$dvs, 2, function(x){(t(x) %*% w0$B %*% x)/sum(abs(s*(D%*%x)))})
    }

    # Generate range of gammas to choose from.
    gammas = sapply(max_gamma, function(x){seq(from=0, to=x, length=num_gammas)})

    ##################################################################################
    # Initialize the validation scores.
    ##################################################################################
    #val_scores = rep(0, times=num_gammas)
    mc_ind = 1
    l0_ind = 1
    best_ind = 1
    min_mc = 1
    min_l0 = p+1

    triv=FALSE

    ##################################################################################
    # Save initial matrices.
    ##################################################################################

    # Initalize objective matrix
    if (dim(w0$B)[2]==1){
      B0 = w0$B
    }    else{
      B0 = t(w0$N) %*% w0$B %*% w0$N
      B0 = (B0+t(B0))/2
    }

    # Initialize the nullspace matrix
    N0 <- w0$N

    # Initialize DVs and iteration lists.
    DVs = list()
    its = list()

    # Initial sparsest solution.
    l0_x = t(N0)%*%t(D)%*%w0$dvs

    # Validation loop
    if(!quiet){
      print("-------------------------------------------")
      print(paste("Fold number:",f))
      print("-------------------------------------------")
    }

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Loop through potential regularization parameters.
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    for(i in 1:num_gammas){
      ##################################################################################
      # Initialization.
      ##################################################################################

      # Initialize output.
      DVs[[i]] = matrix(0, nrow = p, ncol = q)
      its[[i]] = rep(0, times=q)


      # Initialize B and N.
      B = B0
      N = N0

      # Set x0 to be the unpenalized zero-variance discriminant vectors in Null(W0)
      if (dim(B0)[2] == 1){

        # Compute (DN)'*(mu1-mu2)
        w = t(D%*%N0) %*% B0

        # Use normalized w as initial x.
        x0 = w/norm(w,'f')

      }
      else{
        x0 = t(N0)%*% t(D) %*%  w0$dvs[,1]
      }

      # y is the unpenalized solution in the original space.
      # z is the all-zeros vector in the original space.

      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      # Call Alternating Direction Method to solve SZVD problem.
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      for(j in 1:q){
        ## Call ADMM solver.
        tmp = SZVD_ADMM(B = B,  N = N, D=D, pen_scal=s,
                        sols0 = list(x = x0, y = w0$dvs[,j], z= as.matrix(rep(0,p))),
                        gamma=gammas[i,j], beta=beta, tol=tol,
                        maxits=maxits, quiet=TRUE)

        # Extract i-th discriminant vector.
        DVs[[i]][,j] = matrix(D%*%N%*%tmp$x, nrow=p, ncol=1)
        #DVs[[i]][,j] = matrix(tmp$y, nrow=p, ncol=1)

        # Record number of iterations to convergence.
        its[[i]][j] = tmp$its

        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Update N and B for the newly found DV.
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        if(j < (w0$k-1)){
          # Project columns of N onto orthogonal complement of Nx

          x <- as.matrix(DVs[[i]][,j])
          x <- x/norm(as.matrix(x), type = "f")

          # Project N into orthogonal complement of span(x)
          Ntmp <- N - x%*%(crossprod(x,N))

          # Call QR factorization to extract orthonormal basis for span(Ntmp)
          QR = qr(x=Ntmp, LAPACK = TRUE)

          # Extract nonzero rows of R to get columns of Q to use as new N.
          R_rows = (abs(diag(qr.R(QR))) > 1e-6)

          # Use nontrivial columns of Q as updated N.
          N = qr.Q(QR)[, R_rows]

          # Update B0 according to the new basis N.
          B = t(N) %*% w0$B %*% N
          B = 0.5*(B+t(B))

          # Update initial solutions in x direction by projecting next unpenalized ZVD vector.
          x0 = t(N)%*% t(D) %*%  w0$dvs[,(j+1)]
        }

        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Get performance scores on the validation set.
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Call test_ZVD to get predictions, etc.

        statsObj <- test_ZVD(DVs[[i]], Av, w0$means, w0$mu, TRUE)
        stats <- statsObj$stats

        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Validation scores.
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        # Round small entries to zero
        DVs[[i]] <- DVs[[i]] * (ceiling(abs(DVs[[i]])-ztol))

        # if fraction nonzero features less than feat.
        if( 1 <= sum(DVs[[i]] != 0) & sum(DVs[[i]] != 0) <= q*p*feat){
          # Use misclassification rate as validation score.
          scores[f,i] <- stats$mc
        } else if(sum(DVs[[i]] != 0) < 0.5){
          # Trivial solution
          scores[f,i] <- q*p # Disqualify with largest possible value
          triv = TRUE
        } else{
          # Solution is not sparse enough, use most sparse as measure of quality instead.
          scores[f,i] <- sum(DVs[[i]] != 0)
        }

        # Display iteration stats.
        if(!quiet){
          print(sprintf("it = %g, val_score= %g, mc=%g, l0=%g, its=%g", i, scores[f,i],
                        statsObj$stats$mc, sum(statsObj$stats$l0), mean(its[[i]])), quote=F)
        }
      }
      if(triv==TRUE){
        break
      }
    } # End i over validation parameters
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

  ###
  # Solve with lambda = lam(lbest)
  ###
  print(paste("Finished Training: best gamma ind =", gbest))

  # Use the full training set to obtain parameters
  Xt <- X[1:(n-pad),]
  Yt <- Y[1:(n-pad),]
  Atrain <- cbind(Yt%*%matrix(1:K,K,1),Xt)

  ###
  # Generate the gammas
  ###

  # Call ZVD function to solve the unpenalized problem
  w0 <- ZVD(Atrain, scaling = FALSE, get_DVs = TRUE)

  # Extract scaling vector for weighted l1 penalty and diagonal penalty matrix.
  s = sqrt(diag(w0$W))
  w0$s = s

  ## Normalize B (divide by the spectral norm)
  if (dim(w0$B)[2]==1){
    w0$B = w0$B/norm(w0$B, type='f')
  }  else{
    w0$B = (w0$B + t(w0$B))/eigen((w0$B + t(w0$B)), symmetric=TRUE, only.values=TRUE)$values[1]
  }

  # Compute ratio of max gen eigenvalue and l1 norm of the first ZVD to get "bound" on gamma.
  if (dim(w0$B)[2]==1){
    max_gamma =  (t(w0$dvs)%*%w0$B)^2/sum(abs(s*(D %*% w0$dvs)))
  }else{
    max_gamma = apply(w0$dvs, 2, function(x){(t(x) %*% w0$B %*% x)/sum(abs(s*(D%*%x)))})
  }

  # Generate range of gammas to choose from.
  gams = sapply(max_gamma, function(x){seq(from=0, to=x, length=num_gammas)})
  gambest = gammas[gbest,]
  ###

  # Loop until nontrivial solution is found
  trivsol <- TRUE
  while(trivsol){
    szvdObj <- SZVD(Atrain, gambest, D, FALSE, FALSE, tol, maxits, beta, TRUE)
    DVs <- szvdObj$DVs

    # Round small entried to zero
    DVs <- DVs*(ceiling(abs(DVs)-ztol))

    # Check for trivial solution
    if(sum(DVs != 0) == 0){
      # If trivial solution, update gbest by one and update gambest.
      gbest <- gbest - 1
      gambest <- gammas[gbest,]
    } else{
      trivsol <- FALSE
    }
  }

  # Create an object of class SDAPcv to return, might add more to it later
  retOb <- structure(
    list(call = match.call(),
         DVs = DVs,
         gambest = gambest,
         bestInd = gbest),
    class = "SZVD_kFold_cv")

  return(retOb)
}

#' @export
SZVD_kFold_cv.matrix <- function(X, ...){
  res <- SZVD_kFold_cv.default(X, ...)
  cl <- match.call()
  cl[[1L]] <- as.name("SZVD_kFold_cv")
  res$call <- cl
  res
}
