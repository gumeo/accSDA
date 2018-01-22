SDAAPcv <- function (x, ...) UseMethod("SDAAPcv")

SDAAPcv.default <- function(X, Y, folds, Om, gam, lams, q, PGsteps, PGtol, maxits, tol, feat, quiet, initTheta, bt=FALSE, L, eta){
  #
  # HERE WE NEED A DESCRIPTION
  # Use Roxygen2 to create the desired documentation, internal function
  #
  # TODO: handle Y as a factor an generate dummy matrix
  # Get dimensions of input matrices
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

    # Structure to store elements for matrix A used
    # later on, i.e. precomputed values for speed.
    A <- structure(list(
      flag = NA,
      gom = NA,
      X = NA,
      n = NA,
      A = NA
    ),
      class = "Amat"
    )

    # Check if Omega is diagonal
    if(norm(diag(diag(Om))-Om, type = "F") < 1e-15){
      A$flag <- 1
      A$gom <- gam*diag(Om)
      A$X <- Xt
      A$n <- nt
      A$A <- 2*(t(Xt)%*%Xt/nt + gam*Om)
      alpha <- 1/(2*(norm(Xt, type="1")*norm(Xt, type="I")/nt + norm(diag(A$gom), type="I")))
    }else{
      A$flag <- 0
      A$A <- 2*(t(Xt)%*%Xt/nt + gam*Om)
      alpha <- 1/(norm(A$A, type="F"))
    }
    L <- 1/alpha
    L <- norm(diag(diag(Om*gam)),'I')+norm(Xt,'F')^2
    origL <- L
    D <- (1/n)*(t(Yt)%*%Yt)
    R <- chol(D)

    ###
    # Validation loop
    ###
    if(quiet == FALSE){
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
        L <- origL
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

        # Initialize beta
        beta <- matrix(0,p,1)

        ###
        # Alternating direction method to update (theta,beta)
        ###
        for(its in 1:maxits){
          # Compute coefficient vector for elastic net step
          d <- 2*t(Xt)%*%(Yt%*%(theta/nt))

          # Update beta using proximal gradient step
          b_old <- beta
          if(bt == FALSE){
            betaOb <- APG_EN2(A, d, beta, lams[ll], alpha, PGsteps, PGtol)
            beta <- betaOb$x
          }else{
            betaOb <- APG_EN2bt(A, d, beta, lams[ll], L, eta, PGsteps, PGtol)
            #L <- betaOb$L
            beta <- betaOb$x
          }
          # Update theta using projected solution
          if(norm(beta, type = "2") > 1e-12){
            b <- t(Yt)%*%(Xt%*%beta)
            y <- solve(t(R),b)
            z <- solve(R,y)
            tt <- Mj(z)
            t_old <- theta
            theta <- tt/as.numeric(sqrt(t(tt)%*%D%*%tt))

            # Progress
            db <- norm(beta-b_old)/norm(beta, type="2")
            dt <- norm(theta-t_old)/norm(theta, type="2")
          } else{
            # Update beta and theta
            beta <- beta*0
            theta <- theta*0
            # Update changes
            db <- 0
            dt <- 0
          }


          # Check convergence
          if(max(db,dt)<tol){
            # Converged
            break
          }
        }
        # Make the first argument be positive, this is to make the results
        # more reproducible and consistent.
        if(theta[1] < 0){
          theta <- (-1)*theta
          beta <- (-1)*beta
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
  resBest <- SDAAP(Xt, Yt, Om, gam, lams[lbest], q, PGsteps, PGtol, maxits, tol)

  # Create an object of class SDAPcv to return, might add more to it later
  retOb <- structure(
    list(call = match.call(),
         B = resBest$B,
         Q = resBest$Q,
         lbest = lbest,
         lambest = lambest),
    class = "SDAAPcv")

  return(retOb)
}
