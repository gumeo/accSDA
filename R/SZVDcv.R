 #' Cross-validation of sparse zero variance discriminant analysis
#'
#' Applies alternating direction methods of multipliers to solve sparse
#' zero variance descriminant analysis.
#'
#' @param Atrain Training data set.
#' @param Aval Validation set.
#' @param k Number of classes within training and validation sets.
#' @param num_gammas Number of gammas to train on.
#' @param g_mults Parameters defining range of gammas to train, g_max*(c_min, c_max).
#'        Note that it is an array/vector with two elements.
#' @param D Penalty dictionary basis matrix.
#' @param sparsity_pen weight defining validation criteria as weighted sum of misclassification error and
#'        cardinality of discriminant vectors.
#' @param scaling Whether to rescale data so each feature has variance 1.
#' @param penalty Controls whether to apply reweighting of l1-penalty (using sigma = within-class std devs)
#' @param beta Parameter for augmented Lagrangian term in the ADMM algorithm.
#' @param tol Stopping tolerances for the ADMM algorithm, must have tol$rel and tol$abs.
#' @param maxits Maximum number of iterations used in the ADMM algorithm.
#' @param quiet Controls display of intermediate results.
#' @return \code{SZVDcv} returns an object of \code{\link{class}} "\code{SZVDcv}"
#'        including a list with the following named components:
#' \describe{
#'   \item{\code{DVs}}{Discriminant vectors for the best choice of gamma.}
#'   \item{\code{all_DVs}}{Dicriminant vectors for all choices of gamma.}
#'   \item{\code{l0_DVs}}{Discriminant vectors for gamma minimizing cardinality.}
#'   \item{\code{mc_DVs}}{Discriminant vector minimizing misclassification.}
#'   \item{\code{gamma}}{Choice of gamma minimizing validation criterion.}
#'   \item{\code{gammas}}{Set of all gammas trained on.}
#'   \item{\code{max_g}}{Maximum value of gamma guaranteed to yield a nontrivial solution.}
#'   \item{\code{ind}}{Index of best gamma.}
#'   \item{\code{w0}}{unpenalized zero-variance discriminants (initial solutions) plus B and W, etc. from ZVD}
#' }
#' @seealso Non CV version: \code{\link{SZVD}}.
#' @details
#' This function might require a wrapper similar to ASDA.
#' @rdname SZVDcv
#' @export SZVDcv
SZVDcv <- function(Atrain, ...) UseMethod("SZVDcv",Atrain)

#' @return \code{NULL}
#'
#' @rdname SZVDcv
#' @method SZVDcv default
SZVDcv.default <- function(Atrain, Aval, k, num_gammas, g_mults, D, sparsity_pen, scaling = FALSE, penalty = FALSE, beta, tol, maxits, quiet){
  # get dimensions of the training set
  N = dim(Atrain)[1]
  p = dim(Atrain)[2]-1

  ##################################################################################
  ## Compute penalty term for estimating range of regularization parameter values.
  ##################################################################################

  # Call ZVD function to solve the unpenalized problem.
  w0 = ZVD(Atrain, scaling=scaling, get_DVs=TRUE)

  # Extract scaling vector for weighted l1 penalty and diagonal penalty matrix.
  if (penalty==TRUE){ # scaling vector is the std deviations of each feature.
    s = sqrt(diag(w0$W))
  }  else if(penalty==FALSE){ # scaling vector is all-ones (i.e., no scaling)
    s = rep(1, times=p)
  }
  w0$s = s

  # If dictionary D missing, use the identity matri.x
  if (missing(D)){
    D = diag(p)
  }


  ##################################################################################
  ## Compute range of sensible parameter values.
  ##################################################################################

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
  gammas = sapply(max_gamma, function(x){seq(from=g_mults[1]*x, to=g_mults[2]*x, length=num_gammas)})



  ##################################################################################
  ## Get the ZVDs for each choice of gamma and evaluate validation error.
  ##################################################################################

  ##################################################################################
  # Initialize the validation scores.
  ##################################################################################
  val_scores = rep(0, times=num_gammas)
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

  # Initialize nullspace matrix.
  N0 = w0$N

  # Initialize DVs and iteration lists.
  DVs = list()
  its = list()

  # Initial sparsest solution.
  l0_x = t(N0)%*%t(D)%*%w0$dvs


  # y is the unpenalized solution in the original space.
  # z is the all-zeros vector in the original space.

  ##################################################################################
  # For each gamma, calculate ZVDs and corresponding validation scores.

  for (i in (1:num_gammas)){


    ##################################################################################
    # Initialization.
    ##################################################################################

    # Initialize output.
    DVs[[i]] = matrix(0, nrow = p, ncol = (k-1))
    its[[i]] = rep(0, times=(k-1))


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

    ##################################################################################
    ### Get DVs
    ##################################################################################

    for (j in 1:(k-1)){

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


      ##################################################################################
      # Update N and B for the newly found DV.
      ##################################################################################

      if (j < (w0$k-1)) {

        # Project columns of N onto orthogonal complement of Nx.
        x = as.matrix(DVs[[i]][,j])
        x = x/norm(as.matrix(x), 'f')

        # Project N into orthogonal complement of span(x)
        Ntmp = N - x %*% (t(x) %*% N)

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

    }

    ##################################################################################
    # Get performance scores on the validation set.
    ##################################################################################
    # Call test_ZVD to get predictions, etc.
    SZVD_res = test_ZVD(DVs[[i]], Aval, w0$means, w0$mu, scaling=scaling)

    ## Update the cross-validation score for this choice of gamma.

    # If gamma induces the trivial solution, disqualify gamma by assigning
    # large enough penalty that it can't possibly be chosen.
    if (sum(SZVD_res$stats$l0) < 3){
      val_scores[i] = 100*dim(Aval)[1]
      triv=TRUE
    }
    else{
      # if all discriminant vectors are nontrivial use the formula:
      # e_q = %misclassified + % nonzero.
      val_scores[i] = SZVD_res$stats$mc + sparsity_pen*sum(SZVD_res$stats$l0)/(p*(k-1))
    }


    ## Update the best gamma so far.
    # Compare to existing proposed gammas and save best so far.
    if (val_scores[i] <= val_scores[best_ind]){
      best_ind = i
    }

    # Record sparsest nontrivial solution so far.
    if (min(SZVD_res$stats$l0) > 3 && SZVD_res$stats$l0 < min_l0){
      l0_ind = i
      l0_x = DVs[[i]]
      min_l0 = SZVD_res$stats$l0
    }

    # Record best (in terms of misclassification error) so far.
    if (SZVD_res$stats$mc <= min_mc){
      mc_ind = i
      mc_x = DVs[[i]]
      min_mc = SZVD_res$stats$mc
    }


    # Display current iteration stats.
    if (quiet==FALSE){
      print(sprintf("it = %g, val_score= %g, mc=%g, l0=%g, its=%g", i, val_scores[i],
                    SZVD_res$stats$mc, sum(SZVD_res$stats$l0), mean(its[[i]])), quote=F)
    }

    # Terminate if a trivial solution has been found.
    if (triv==TRUE){
      break
    }

  }
  ##################################################################################

  # Export discriminant vectors found using validation.
  val_x = DVs[[best_ind]]

  # Return best ZVD, gamma, lists of gammas and validation scores, etc.
  return(structure(list(call = match.call(),
                        DVs = val_x,
                        all_DVs = DVs,
                        l0_DVs = l0_x,
                        mc_DVs = mc_x,
                        gamma = gammas[best_ind,],
                        gammas = gammas,
                        max_g = max_gamma,
                        ind=best_ind,
                        scores=val_scores,
                        w0=w0,
                        x0=x0),
              class="SZVDcv"))
}

#' @export
SZVDcv.matrix <- function(Atrain, ...){
  res <- SZVDcv.default(Atrain, ...)
  cl <- match.call()
  cl[[1L]] <- as.name("SZVDcv")
  res$call <- cl
  res
}
