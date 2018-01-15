SDADcv <- function (x, ...) UseMethod("SDADcv")

SDADcv.default <- function(X, Y, folds, Om, gam, lams, mu, q, PGsteps, PGtol, maxits, tol, feat, quiet, initTheta){
  #
  # HERE WE NEED A DESCRIPTION
  # Use Roxygen2 to create the desired documentation
  #
  # TODO: handle Y as a factor an generate dummy matrix
  # Get dimensions of input matrices
  # Note: PGtol is a vector with abs(index 1) and rel(index 2) tolerances
  # Inconsistency with SDAD on norm of beta, here 1e-12, 1e-15 in SDAD
  dimX <- dim(X)
  n <- dimX[1]
  p <- dimX[2]
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

  # Sort lambdas in descending order
  lams <- lams[order(lams,decreasing = TRUE)]

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
  nlam <- length(lams)

  # Validation scores
  scores <- q*p*matrix(1,nrow = folds, ncol = nlam)

  # Misclassification rate for each classifier
  mc <- matrix(0,nrow = folds, ncol = nlam)

  for(f in 1:folds){
    # Initialization

    # Extract X and Y from training data
    Xt <- X[tinds,]
    Yt <- Y[tinds,]

    # Extract validation data
    Xv <- X[vinds,]
    Yv <- Y[vinds,]

    # Get dimensions of training matrices
    nt <- dim(Xt)[1]
    p  <- dim(Xt)[2]

    # Centroid matrix of training data
    C <- diag(diag((1/(t(Yt)%*%Yt))))%*%t(Yt)%*%Xt

    # Check if Om is diagonal. If so, use matrix inversion lemma in linear
    # system solves.
    if(norm(diag(diag(Om))-Om,type = "F") < 1e-15){
      # Flag to use Sherman-Morrison-Woodbury to translate to
      # smaller dimensional linear system solves.
      SMW <- 1

      # Easy to invert diagonal part of Elastic net coefficient matrix.
      M <- mu*diag(p) + 2*gam*Om
      Minv = 1/diag(M)

      # Cholesky factorization for smaller linear system.
      RS = chol(diag(nt) + 2*Xt%*%diag(Minv)%*%t(Xt)/nt);
    } else{ # Use Cholesky for solving linear systems in ADMM step
      # Flag to not use SMW
      SMW <- 0
      A <- mu*diag(p) + 2*(t(Xt)%*%Xt + gam*Om)
      R2 <- chol(A)
    }
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Matrices for theta update.
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    D <- (1/nt)*(t(Yt)%*%Yt)
    R <- chol(D) # Cholesky factorization of D.

    ###
    # Validation loop
    ###
    if(!quiet){
      print("-------------------------------------------")
      print(paste("Fold number:",f))
      print("-------------------------------------------")
    }

    ###
    # Loop through the validation parameters
    ###
    for(ll in 1:nlam){
      # Initialize B and Q
      Q <- matrix(1,K,q)
      B <- matrix(0,p,q)

      #-------------------------------------------------
      # Call Alternating Direction Method to solve SDA
      #-------------------------------------------------
      # For j=1,2,...,q compute the SDA pair (theta_j, beta_j)
      for(j in 1:q){
        # Initialization

        # Compute Qj (K by j, first j-1 scoring vectors, all-ones last col)
        Qj <- Q[,1:j]

        # Precompute Mj = I-Qj*Qj'*D
        Mj <- function(u){
          return(u-Qj%*%(t(Qj)%*%(D%*%u)))
        }

        # Initialize theta
        theta <- matrix(stats::runif(K),nrow=K,ncol=1)
        theta <- Mj(theta)
        if(j == 1 & !missing(initTheta)){
          theta=initTheta
        }
        theta <- theta/as.numeric(sqrt(crossprod(theta,D%*%theta)))

        # Initialize coefficient vector for elastic net step
        d <- 2*t(Xt)%*%(Yt%*%theta)

        # Initialize beta
        if(SMW == 1){
          btmp <- Xt%*%(Minv*d)/nt
          beta <- (Minv*d) - 2*Minv*(t(Xt)%*%(solve(RS,solve(t(RS),btmp))))
        }else{
          beta <- solve(R2,solve(t(R2),d))
        }

        ###
        # Alternating direction method to update (theta,beta)
        ###
        for(its in 1:maxits){
          # Update beta using alternating direction method of multipliers.
          b_old <- beta

          if(SMW == 1){
            # Use SMW-based ADMM
            betaOb <- ADMM_EN_SMW(Minv, Xt, RS, d, beta, lams[ll], mu, PGsteps, PGtol, TRUE)
            beta <- betaOb$y
          } else{
            betaOb <- ADMM_EN2(R2, d, beta, lams[ll], mu, PGsteps, PGtol, TRUE)
            beta <- betaOb$y
          }

          # Update theta using the projected solution
          if(norm(beta, type="2") > 1e-12){
            # Update theta
            b <- t(Yt)%*%(Xt%*%beta)
            y <- solve(t(R),b)
            z <- solve(R,y)
            tt <- Mj(z)
            t_old <- theta
            theta <- tt/sqrt(as.numeric(t(tt)%*%D%*%tt))

            # Update changes
            db <- norm(beta-b_old, type="2")/norm(beta, type="2")
            dt <- norm(theta-t_old, type="2")/norm(theta, type="2")
          } else{
            # Update b and theta
            beta <- beta*0
            theta <- theta*0
            db <- 0
            dt <- 0
          }

          if(max(db,dt) < tol){
            # Converged
            break
          }
        }
        # Update Q and B
        Q[,j] <- theta
        B[,j] <- beta
      }

      #------------------------------------------------------------
      # Get classification statistics for (Q,B)
      #------------------------------------------------------------

      # Project validation data
      PXtest <- Xv%*%B
      # Project centroids
      PC <- C%*%B

      # Compute distances to the centroid for each projected test observation
      dists <- matrix(0,nv,K)
      for(i in 1:nv){
        for(j in 1:K){
          dists[i,j] <- norm(PXtest[i,] - PC[j,], type="2")
        }
      }

      # Label test observation according to the closest centroid to its projection.
      predicted_labels <- t(apply(dists, 1, function(x) c(min(x),which.min(x))))
      predicted_labels <- predicted_labels[,2] # Select the indices

      # Form predicted Y
      Ypred <- matrix(0,nv,K)
      for(i in 1:nv){
        Ypred[i,predicted_labels[i]] <- 1
      }

      # Fraction misclassified
      mc[f,ll] <- (0.5*norm(Yv-Ypred,type="F")^2)/nv

      ###
      # Validation scores
      ###
      # if fraction nonzero features less than feat.
      if( 1 <= sum(B != 0) & sum(B != 0) <= q*p*feat){
        # Use misclassification rate as validation score.
        scores[f,ll] <- mc[f,ll]
      } else if(sum(B != 0) > q*p*feat){
        # Solution is not sparse enough, use most sparse as measure of quality instead.
        scores[f,ll] <- sum(B != 0)
      }

      # Display iteration stats
      if(!quiet){
        print(paste("f:", f, "| ll:", ll, "| lam:", lams[ll], "| feat:",
                    sum(B != 0)/(q*p), "| mc:", mc[f,ll], "| score:", scores[f,ll]))
      }
    } # End of for ll in 1:nlam
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
  } # End of folds loop
  ###
  # Find the best solution
  ###

  # Average CV scores
  avg_score <- colMeans(scores)

  # Choose lambda with best average score
  lbest <- which.min(avg_score)

  lambest <- lams[lbest]

  ###
  # Solve with lambda = lam(lbest)
  ###
  print(paste("Finished Training: lam =", lambest))

  # Use the full training set to obtain parameters
  Xt <- X[1:(n-pad),]
  Yt <- Y[1:(n-pad),]

  # Get best Q and B on full training data
  resBest <- SDAD(Xt, Yt, Om, gam, lambest, mu, q, PGsteps, PGtol, maxits, tol)

  # Create an object of class SDAPcv to return, might add more to it later
  retOb <- structure(
    list(call = match.call(),
         B = resBest$B,
         Q = resBest$Q,
         lbest = lbest,
         lambest = lambest),
    class = "SDADcv")

  return(retOb)

}
