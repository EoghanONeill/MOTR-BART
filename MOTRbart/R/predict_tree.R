#' @export
predict_motr_bart = function(object, newdata,
                             type = c('all', 'median', 'mean')) {
  # Get the means and sds to standardise the covariates from the test data

  center = object$center_x
  scale = object$scale_x
  newdata = as.matrix(cbind(1,scale(newdata, center=center, scale=scale)))
  ancestors = object$ancestors

  # Create holder for predicted values
  n_newX = dim(newdata)[1]
  n_its = object$npost
  y_hat_mat = matrix(NA, nrow = n_its,
                     ncol = nrow(newdata))
  num_tress = object$num_trees

  # Now loop through iterations and get predictions
  for (i in 1:n_its) {
    # Get current set of trees
    curr_trees = object$trees[[i]]

    # Use get_predictions function to get predictions
    y_hat_mat[i,] = get_predictions(curr_trees,
                                    newdata,
                                    single_tree = length(curr_trees) == 2,
                                    ancestors = ancestors)
  }

  # Sort out what to return
  out = switch(type,
               all = object$y_mean + object$y_sd * (y_hat_mat + (object$y_max + object$y_min)/2),
               mean = object$y_mean + object$y_sd * (apply(y_hat_mat,2,'mean')+ (object$y_max + object$y_min)/2),
               median = object$y_mean + object$y_sd * (apply(y_hat_mat,2,'median')+ (object$y_max + object$y_min)/2))

  return(out)

} # end of predict function


#' @export
predict_TVPbart = function(object, newdata,
                             type = c('all', 'median', 'mean')) {
  # Get the means and sds to standardise the covariates from the test data

  center = object$center_x
  scale = object$scale_x
  newdata = as.matrix(cbind(1,scale(newdata, center=center, scale=scale)))
  ancestors = object$ancestors

  # Create holder for predicted values
  n_newX = dim(newdata)[1]
  n_its = object$npost
  y_hat_mat = matrix(NA, nrow = n_its,
                     ncol = nrow(newdata))
  num_tress = object$num_trees

  Ltest <- matrix(1,nrow = nrow(newdata), ncol = object$ntrain)

  # Now loop through iterations and get predictions
  for (i in 1:n_its) {
    # Get current set of trees
    curr_trees = object$trees[[i]]

    # Use get_predictions function to get predictions
    y_hat_mat[i,] = TVPget_predictions(curr_trees, Ltest,
                                    newdata,
                                    single_tree = length(curr_trees) == 2) +
      rnorm(n_newX,mean = 0, sd = sqrt(object$vars_betas[i,][1]))
  }

  # Sort out what to return
  out = switch(type,
               all = object$y_mean + object$y_sd * (y_hat_mat + (object$y_max + object$y_min)/2),
               mean = object$y_mean + object$y_sd * (apply(y_hat_mat,2,'mean')+ (object$y_max + object$y_min)/2),
               median = object$y_mean + object$y_sd * (apply(y_hat_mat,2,'median')+ (object$y_max + object$y_min)/2))

  return(out)

} # end of predict function


########################################################################################################
# Predictions for classification
########################################################################################################

#' @export
predict_motr_bart_class = function(object,
                                   newdata,
                                   type = c('all', 'median', 'mean')) {
  # Get the means and sds to standardise the covariates from the test data
  center = object$center_x
  scale = object$scale_x
  newdata = as.matrix(cbind(1,scale(newdata, center=center, scale=scale)))
  ancestors = object$ancestors

  # Create holder for predicted values
  n_newX = dim(newdata)[1]
  n_its = object$npost
  y_hat_mat = matrix(NA, nrow = n_its,
                     ncol = nrow(newdata))
  num_tress = object$num_trees

  # Now loop through iterations and get predictions
  for (i in 1:n_its) {
    # Get current set of trees
    curr_trees = object$trees[[i]]

    # Use get_predictions function to get predictions
    y_hat_mat[i,] = get_predictions(curr_trees,
                                    newdata,
                                    single_tree = length(curr_trees) == 1,
                                    ancestors = ancestors)
  }

  # Sort out what to return
  out = switch(type,
               all = y_hat_mat,
               mean = apply(pnorm(y_hat_mat),2,'mean'),
               median = apply(pnorm(y_hat_mat),2,'median'))

  return(out)

} # end of predict function

