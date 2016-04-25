#' Classify test data using nearest centroid classification and
#' discriminant vectors leatned from the training set.
#'
#' This function is used in SZVDcv and is only meant for internal
#' use at this stage. Will potentially be released in future versions.
#'
#' @param w Matrix with columns equal to discriminant vectors.
#' @param test matrix containing test set.
#' @param classMeans Means of each class in the training set,
#'        (used for computing centroids for classification).
#' @param mus means/standard devs of the training set,
#'        (used for centering/normalizing the test data appropriately).
#' @param scaling Logical indicating wether scaling should be done.
#'        on the test set.
#' @return \code{test_ZVD} returns an object of \code{\link{class}} "\code{test_ZVD}" including a list
#' with the following named components
#'
#' \describe{
#'   \item{\code{stats}}{list containing number of misclassified observations,
#'   l0 and l1 norms of discriminants.}
#'   \item{\code{pred}}{predicted class labels according to nearest centroid
#'   and the discriminants.}
#' }
#' @seealso Used by: \code{\link{SZVDcv}}.
#' @details
#' This function is used by other functions and should only be called explicitly for
#' debugging purposes. Potential release in the future.
#' @keywords internal
test_ZVD <- function(w, test, classMeans, mus, scaling){
  # 
  if(scaling){
    mu <- as.matrix(mus$mu)
    sig <- mus$sig
  } else{
    mu <- as.matrix(mus)
  }
  
  # Extract class labels and observations
  test_labels <- test[,1]
  test_obs <- test[,2:dim(test)[2]]
  
  # Get number of test observations
  N <- length(test_labels)
  
  # Get the number of classes
  K <- max(test_labels)
  
  # Center the test data
  test_obs <- test_obs - matrix(1,N,1)%*%t(mu)
  
  # Scale according to the saved scaling from the training data (if desired)
  if(scaling){
    test_obs <- test_obs %*% diag(1/sig)
  }
  
  #====================================================================
  # Classify the test data according to nearest centroid rule.
  #====================================================================
    
  # Project the test data to the lower dim linear space defined by the ZVDs.
  proj <- t(w)%*%t(test_obs)
  
  # Compute centroids and projected distances.
  cent <- crossprod(w,classMeans)
  
  # Compute distances to the centroid for each projected test observation.
  dists <- matrix(0,N,K)
  for(i in 1:N){
    for(j in 1:K){
      dists[i,j] <- norm(proj[,i] - cent[,j], type = "2")
    }
  }
  
  # Label test observation according to the closes centroid in its projection.
  predicted_labels <- t(apply(dists, 1, function(x) c(min(x),which.min(x))))
  predicted_labels <- predicted_labels[,2] # Select the indices
  
  #===================================================================
  # Compute fraction misclassed, l0 and l1 norms of the classifiers.
  #====================================================================
    
  # Compute fraction of misclassified observations.
  misclassed <- sum(abs(test_labels - predicted_labels) > 0) / N
  
  # l0
  l0 <- sum(abs(w) > 1e-3)
  
  # l1
  l1 <- sum(abs(w))
  
  #====================================================================
  # Output results.
  #====================================================================
  
  # stats
  stats <- list(mc = misclassed,
                l0 = l0,
                l1 = l1)
  
  # Create an object of class test_ZVD to return
  retOb <- structure(
    list(call = match.call(),
         stats = stats,
         preds = predicted_labels),
    class = "test_ZVD")
  
  return(retOb)
}