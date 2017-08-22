#' Classify test data using nearest centroid classification and
#' discriminant vectors learned from the training set.
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
#' @param scaling Logical indicating whether scaling should be done.
#'        on the test set.
#' @param ztol Threshold for setting values in DVs to zero.
#' @return \code{test_ZVD} returns an object of \code{\link{class}} "\code{test_ZVD}" including a list
#' with the following named components
#'
#' \describe{
#'   \item{\code{stats}}{list containing number of misclassified observations,
#'   l0 and l1 norms of discriminants.}
#'   \item{\code{pred_labs}}{predicted class labels according to nearest centroid
#'   and the discriminants.}
#' }
#' @seealso Used by: \code{\link{SZVDcv}}.
#' @details
#' This function is used by other functions and should only be called explicitly for
#' debugging purposes. Potential release in the future.
#' This function should potentially be made internal for the release.
#' @rdname test_ZVD
#' @export test_ZVD
test_ZVD <- function(w, test, classMeans, mus, scaling, ztol){
  ######################################################################################
  # Initialization.
  ######################################################################################

  # Get scaling/centering factors.
  if (scaling==TRUE){
    mu = mus$mu
    sig = mus$sig
  } else{
    mu = mus
  }


  # Extract class labels and observations.
  test_labels = factor(test[,1])
  test_obs = as.matrix(data.frame(test[,2:dim(test)[2]]))

  # Get number of test observations.
  N = length(test_labels)

  # Get number of classes.
  K = length(levels(test_labels))

  # Center the test data.
  test_obs = test_obs - matrix(1, nrow =N, ncol=1) %*%  t(mu)

  # Scale according to the saved scaling from the training data (if desired)
  if (scaling==TRUE){
    test_obs = test_obs %*% diag(1/sig)
  }

  ######################################################################################
  # Classify the test data according to nearest centroids.
  ######################################################################################


  ## Project the test data to the lower dim linear space defined by the ZVDs.
  proj = t(w) %*% t(as.matrix(test_obs))


  ## Compute centroids and projected distances.
  cent = t(w) %*% classMeans

  # Compute distances to the centroid for test data.
  dist = apply(X = t(proj), MARGIN=1,
               FUN=function(y){
                 apply(X=cent, MARGIN=2, FUN= function(x){ norm(as.matrix( x - y, 'f')) } )
               }
  )

  # Label test observation according to the closest centroid to its projection.
  predicted_labels = max.col(-t(dist))

  # Extract labels from the test data (assumed these are integers from 1,..., k.)
  true_labels = test[,1]



  ######################################################################################
  # Compute misclassed, l0 and l1 norms of the classifiers.
  ######################################################################################

  # Compute fraction of misclassified observations.
  misclassed= sum(abs(true_labels - predicted_labels) > 0) / N

  # l0
  l0 = apply(w, MARGIN=2, FUN= function(x){ sum(abs(x)>ztol)})

  # l1
  l1 = apply(w, MARGIN=2, FUN= function(x){sum(abs(x))})


  ######################################################################################
  # Output results.
  ######################################################################################
  results = list(stats=list(mc=misclassed, l0=l0, l1=l1), pred_labs = predicted_labels, dist=dist)
  return(results)
}
