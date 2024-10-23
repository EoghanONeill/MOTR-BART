
#' @import collapse
#' @import SoftBart
#' @import mvnfast
#' @importFrom mvtnorm 'rmvnorm'
#' @importFrom stats 'rgamma' 'runif' 'dnorm' 'sd' 'rnorm' 'pnorm'
#' @importFrom MCMCpack 'rdirichlet'
#' @importFrom truncnorm 'rtruncnorm'
#' @importFrom Rfast 'sort_unique'
#' @export

motr_bart_alpha = function(x,
                           y,
                           sparse = TRUE,
                           alpha_a = 0.5, # Linero alpha prior parameter
                           alpha_b = 1, # Linero alpha prior parameter
                           vars_inter_slope = TRUE,
                           ntrees = 10,
                           node_min_size = 5,
                           alpha = 0.95,
                           beta = 2,
                           nu = 3,
                           lambda = 0.1,
                           sigma2 = 1,
                           nburn = 1000,
                           npost = 1000,
                           nthin = 1,
                           ancestors = FALSE,
                           trans_prob = c(2.5, 2.5, 4) / 9, # Probabilities to grow, prune or change, respectively
                           alpha_prior = FALSE,
                           max_bad_trees = 10,
                           splitting_rules = "discrete",
                           coeff_prior_conj = TRUE,
                           k = 2,
                           sigquant = .90,
                           centre_y = TRUE) {


  if(nrow(x) != length(y)){
    stop("nrow(x) != length(y)")
  }

  if(!(splitting_rules %in% c("discrete", "continuous"))){
    stop("splitting_rules must be 'discrete' or 'continuous'.")
  }


  x = as.data.frame(x)
  # Quantities needed for prediction
  center = apply(x, 2, mean)
  scale = apply(x, 2, sd)
  aux.X = lapply(x, unique) # Checking how many unique values each variable has
  unique.values.X = unlist(lapply(aux.X, length))

  center[which(unique.values.X<=2)] = 0
  scale[which(unique.values.X<=2)] = 1

  X_orig = x
  X = as.matrix(cbind(1,scale(x, center, scale))) # standardising the covariates and adding an intercept


  # print("line 45 ncol(X) = ")
  # print(ncol(X))
  #
  # print("nrow(X) = ")
  # print(nrow(X))
  # Extract control parameters
  node_min_size = node_min_size

  # Extract MCMC details
  TotIter = nburn + npost*nthin # Total of iterations

  # Storage containers
  store_size = npost
  tree_store = vector('list', store_size)
  sigma2_store = rep(NA, store_size)
  y_hat_store = matrix(NA, ncol = length(y), nrow = store_size)
  var_count = rep(0, ncol(X_orig))
  var_count_store = matrix(0, ncol = ncol(X_orig), nrow = store_size)
  s_prob_store = matrix(0, ncol = ncol(X_orig), nrow = store_size)
  vars_betas_store = matrix(0, ncol = 2, nrow = store_size)
  tree_fits_store = matrix(0, ncol = ntrees, nrow = length(y))
  alpha_s_store <- rep(NA, store_size)

  # Scale the response target variable
  y_mean = mean(y)
  y_sd = sd(y)
  y_scale = (y - y_mean)/y_sd
  n = length(y_scale)
  p = ncol(X_orig)
  s = rep(1/p, p)

  # Prior for the beta vector

  if(centre_y){
    y_max <- max(y_scale)
    y_min <- min(y_scale)
  }else{
    y_max <- 0
    y_min <- 0
  }
  y_scale <- y_scale - (y_max + y_min)/2

  # sigma2_beta <- (max(y_scale)-min(y_scale))/((2 * k * sqrt(ntrees))^2)
  sigma2_beta <- ((max(y_scale)-min(y_scale))/(2 * k * sqrt(ntrees)))^2

  tau_b = ntrees

  if(coeff_prior_conj == FALSE){
    tau_b <- 1/sigma2_beta
  }
  V = rep(1/tau_b, 2)
  inv_V = 1/V


  if(alpha_prior == TRUE){
    alpha_s <- p
  }else{
    alpha_s <- 1
  }
  alpha_scale <- p








  if( coeff_prior_conj == FALSE ) {
    if(p < n) {
      df = data.frame(x,y_scale)
      lmf = lm(y_scale~.,df)
      sigest = summary(lmf)$sigma
    } else {
      sigest = sd(y_scale)
    }

    qchi = qchisq(1.0-sigquant,nu)
    lambda = (sigest*sigest*qchi)/nu #lambda parameter for sigma prior


    sigma2 <- sigest^2
  }

  # print("line 128")

  # s = rep(1/p,p)


  # Create a list of trees for the initial stump
  curr_trees = create_stump(num_trees = ntrees,
                            y = y_scale,
                            X = X)
  # Initialise the new trees as current one
  new_trees = curr_trees

  # Initialise the predicted values to zero
  predictions = get_predictions(curr_trees, X, single_tree = ntrees == 1, ancestors)

  # Set up a progress bar
  pb = utils::txtProgressBar(min = 1, max = TotIter,
                             style = 3, width = 60,
                             title = 'Running rBART...')

  # Start the MCMC iterations loop
  for (i in 1:TotIter) {

    utils::setTxtProgressBar(pb, i)



    # Start looping through trees
    for (j in 1:ntrees) {

      current_partial_residuals = y_scale - predictions + tree_fits_store[,j]

      # Propose a new tree via grow/change/prune/swap
      # type = sample(c('grow', 'prune', 'change'), 1)

      type = sample_move(curr_trees[[j]], i, nburn,
                         trans_prob)

      # if(get_nterminal(curr_trees[[j]]) ==1){ type = 'grow'} # Grow the tree if it is a stump
      # if(i < max(floor(0.1*nburn), 5)){ type = 'grow'} # Grow for the first few iterations

      # print("ncol(X) = ")
      # print(ncol(X))
      #
      # print("nrow(X) = ")
      # print(nrow(X))

      # Generate a new tree based on the current
      new_trees[[j]] = update_tree(y = y_scale,
                                   X = X,
                                   type = type,
                                   curr_tree = curr_trees[[j]],
                                   node_min_size = node_min_size,
                                   s = s,
                                   max_bad_trees = max_bad_trees,
                                   splitting_rules = splitting_rules)

      # CURRENT TREE: compute the log of the marginalised likelihood + log of the tree prior
      l_old = tree_full_conditional(curr_trees[[j]],
                                    X,
                                    current_partial_residuals,
                                    sigma2,
                                    V,
                                    inv_V,
                                    nu,
                                    lambda,
                                    tau_b,
                                    ancestors,
                                    coeff_prior_conj) # + get_tree_prior(curr_trees[[j]], alpha, beta)

      # NEW TREE: compute the log of the marginalised likelihood + log of the tree prior
      l_new = tree_full_conditional(new_trees[[j]],
                                    X,
                                    current_partial_residuals,
                                    sigma2,
                                    V,
                                    inv_V,
                                    nu,
                                    lambda,
                                    tau_b,
                                    ancestors,
                                    coeff_prior_conj) # + get_tree_prior(new_trees[[j]], alpha, beta)

      # # Exponentiate the results above
      # if(type == 'grow'){
      #   # a = exp(l_new - l_old)*ratio_grow(new_trees[[j]], curr_trees[[j]])
      #   a = exp(l_new - l_old)*ratio_grow(curr_trees[[j]], new_trees[[j]])
      # } else if(type == 'prune'){
      #   # a = exp(l_new - l_old)*ratio_prune(new_trees[[j]], curr_trees[[j]])
      #   a = exp(l_new - l_old)*ratio_prune(curr_trees[[j]], new_trees[[j]])
      # } else{
      #   a = exp(l_new - l_old)
      # }


      a <- get_MH_probability2(curr_trees[[j]], new_trees[[j]], l_old, l_new,
                               type, trans_prob,
                               alpha, beta)




      if(a > runif(1)) {

        curr_trees[[j]] = new_trees[[j]] # The current tree "becomes" the new tree, if the latter is better

        if (type =='change'){
          var_count[curr_trees[[j]]$var[1] - 1] = var_count[curr_trees[[j]]$var[1] - 1] - 1
          var_count[curr_trees[[j]]$var[2] - 1] = var_count[curr_trees[[j]]$var[2] - 1] + 1
        }

        if (type=='grow'){
          var_count[curr_trees[[j]]$var - 1] = var_count[curr_trees[[j]]$var - 1] + 1 } # -1 because of the intercept in X

        if (type=='prune'){
          var_count[curr_trees[[j]]$var - 1] = var_count[curr_trees[[j]]$var - 1] - 1 } # -1 because of the intercept in X
      }

      # Update mu whether tree accepted or not
      curr_trees[[j]] = simulate_beta(curr_trees[[j]],
                                      X,
                                      current_partial_residuals,
                                      sigma2,
                                      inv_V,
                                      tau_b,
                                      nu,
                                      ancestors,
                                      coeff_prior_conj)

      current_fit = get_predictions(curr_trees[[j]], X, single_tree = TRUE, ancestors)
      predictions = predictions - tree_fits_store[,j] # subtract the old fit
      predictions = predictions + current_fit # add the new fit
      tree_fits_store[,j] = current_fit # update the new fit

    } # End loop through trees

    # Updating the predictions (y_hat)
    # predictions = get_predictions(curr_trees, X, single_tree = ntrees == 1, ancestors)

    sum_of_squares = sum((y_scale - predictions)^2)

    # Update sigma2 (variance of the residuals)
    sigma2 = update_sigma2(sum_of_squares, n = length(y_scale), nu, lambda)

    # Update sigma2_beta0 and sigma2_beta1
    if (vars_inter_slope == TRUE) {
      vars_betas = update_vars_intercepts_slopes(curr_trees, ntrees, sigma2, coeff_prior_conj)
      V = 1/c(vars_betas$var_inter, vars_betas$var_slopes)
      inv_V = 1/V
    }

    # # Update s = (s_1, ..., s_p), where s_p is the probability that predictor p is used to create new terminal nodes
    # if (sparse == 'TRUE' & i > floor(TotIter*0.1)){
    #   s = update_s(var_count, p, 1)
    # }


    # Update s = (s_1, ..., s_p), where s_p is the probability that predictor q in 1:p is used to create new terminal nodes
    if (sparse & (i > floor(nburn * 0.5))) {
      s_update <- update_s(var_count, p, alpha_s)
      s <- s_update[[1]]
      if(length(s) != ncol(X_orig)){
        print("s = ")
        print(s)
        print("length(s) = ")
        print(length(s))
        print("ncol(X) = ")
        print(ncol(X_orig))
        print("p = ")
        print(p)

      }

      if(alpha_prior == TRUE){
        alpha_s <- update_alpha(s, alpha_scale, alpha_a, alpha_b, p, s_update[[2]])
      }
    }

    # If at the right place, store everything
    if((i > nburn) & ((i - nburn) %% nthin) == 0) {
      curr = (i - nburn)/nthin
      tree_store[[curr]] = curr_trees
      sigma2_store[curr] = sigma2
      y_hat_store[curr,] = predictions
      var_count_store[curr,] = var_count
      s_prob_store[curr,] = s
      vars_betas_store[curr,] = V
      alpha_s_store[curr] <- alpha_s
    }

  } # End iterations loop

  cat('\n') # Make sure progress bar ends on a new line

  return(list(trees = tree_store,
              sigma2 = sigma2_store*y_sd^2,
              y_hat = (y_hat_store + (y_max + y_min)/2 )*y_sd + y_mean,
              center_x = center,
              scale_x = scale,
              npost = npost,
              nburn = nburn,
              nthin = nthin,
              ntrees = ntrees,
              y_mean = y_mean,
              y_sd = y_sd,
              ancestors = ancestors,
              var_count_store = var_count_store,
              s = s_prob_store,
              vars_betas = vars_betas_store,
              alpha_s_store = alpha_s_store,
              y_max = y_max,
              y_min = y_min))

} # End main function




#' @import collapse
#' @import SoftBart
#' @import mvnfast
#' @importFrom mvtnorm 'rmvnorm'
#' @importFrom stats 'rgamma' 'runif' 'dnorm' 'sd' 'rnorm' 'pnorm'
#' @importFrom MCMCpack 'rdirichlet'
#' @importFrom truncnorm 'rtruncnorm'
#' @importFrom Rfast 'sort_unique'
#' @export

TVPbart = function(x,
                           y,
                           sparse = TRUE,
                           alpha_a = 0.5, # Linero alpha prior parameter
                           alpha_b = 1, # Linero alpha prior parameter
                           vars_inter_slope = TRUE,
                           ntrees = 10,
                           node_min_size = 5,
                           alpha = 0.95,
                           beta = 2,
                           nu = 3,
                           lambda = 0.1,
                           sigma2 = 1,
                           nburn = 1000,
                           npost = 1000,
                           nthin = 1,
                           ancestors = FALSE,
                           trans_prob = c(2.5, 2.5, 4) / 9, # Probabilities to grow, prune or change, respectively
                           alpha_prior = FALSE,
                           max_bad_trees = 10,
                           splitting_rules = "discrete",
                   coeff_prior_conj = TRUE,
                   k = 2,
                   sigquant = .90,
                   centre_y = TRUE) {


  if(nrow(x) != length(y)){
    stop("nrow(x) != length(y)")
  }

  if(!(splitting_rules %in% c("discrete", "continuous"))){
    stop("splitting_rules must be 'discrete' or 'continuous'.")
  }


  x = as.data.frame(x)
  # Quantities needed for prediction
  center = apply(x, 2, mean)
  scale = apply(x, 2, sd)
  aux.X = lapply(x, unique) # Checking how many unique values each variable has
  unique.values.X = unlist(lapply(aux.X, length))

  center[which(unique.values.X<=2)] = 0
  scale[which(unique.values.X<=2)] = 1

  X_orig = x
  X = as.matrix(cbind(1,scale(x, center, scale))) # standardising the covariates and adding an intercept


  # print("line 45 ncol(X) = ")
  # print(ncol(X))
  #
  # print("nrow(X) = ")
  # print(nrow(X))
  # Extract control parameters
  # node_min_size = node_min_size

  # Extract MCMC details
  TotIter = nburn + npost*nthin # Total of iterations

  # Storage containers
  store_size = npost
  tree_store = vector('list', store_size)
  sigma2_store = rep(NA, store_size)
  y_hat_store = matrix(NA, ncol = length(y), nrow = store_size)
  var_count = rep(0, ncol(X_orig))
  var_count_store = matrix(0, ncol = ncol(X_orig), nrow = store_size)
  s_prob_store = matrix(0, ncol = ncol(X_orig), nrow = store_size)
  vars_betas_store = matrix(0, ncol = 1, nrow = store_size)
  tree_fits_store = matrix(0, ncol = ntrees, nrow = length(y))
  alpha_s_store <- rep(NA, store_size)

  # Scale the response target variable
  y_mean = mean(y)
  y_sd = sd(y)
  y_scale = (y - y_mean)/y_sd
  n = length(y_scale)
  p = ncol(X_orig)
  s = rep(1/p, p)

  if(centre_y){
    y_max <- max(y_scale)
    y_min <- min(y_scale)
  }else{
    y_max <- 0
    y_min <- 0
  }
  y_scale <- y_scale - (y_max + y_min)/2

  # Prior for the beta vector
  # sigma2_beta <- (max(y_scale)-min(y_scale))/((2 * k * sqrt(ntrees))^2)

  sigma2_beta <- (max(y_scale)-min(y_scale))/((2 * k * sqrt(ntrees*p))^2)

  tau_b = 1 #ntrees*p

  if(coeff_prior_conj == FALSE){
    tau_b <- 1/sigma2_beta
  }

  V = rep(1/tau_b, 1)
  inv_V = 1/V


  if(alpha_prior == TRUE){
    alpha_s <- p
  }else{
    alpha_s <- 1
  }
  alpha_scale <- p


  # print("line 128")

  # s = rep(1/p,p)
  Lmatleaf <- 1*lower.tri(diag(n),diag = TRUE)


  # Create a list of trees for the initial stump
  curr_trees = TVPcreate_stump(num_trees = ntrees,
                            y = y_scale,
                            Lmat = Lmatleaf)
  # Initialise the new trees as current one
  new_trees = curr_trees

  # Initialise the predicted values to zero
  # predictions = TVPget_predictions(curr_trees, Lmatleaf, X, single_tree = ntrees == 1)

  predictions = numeric(n)




  # Set up a progress bar
  pb = utils::txtProgressBar(min = 1, max = TotIter,
                             style = 3, width = 60,
                             title = 'Running rBART...')

  # Start the MCMC iterations loop
  for (i in 1:TotIter) {

    utils::setTxtProgressBar(pb, i)



    # Start looping through trees
    for (j in 1:ntrees) {

      current_partial_residuals = y_scale - predictions + tree_fits_store[,j]

      # Propose a new tree via grow/change/prune/swap
      # type = sample(c('grow', 'prune', 'change'), 1)

      type = sample_move(curr_trees[[j]], i, nburn,
                         trans_prob)

      # if(get_nterminal(curr_trees[[j]]) ==1){ type = 'grow'} # Grow the tree if it is a stump
      # if(i < max(floor(0.1*nburn), 5)){ type = 'grow'} # Grow for the first few iterations

      # print("ncol(X) = ")
      # print(ncol(X))
      #
      # print("nrow(X) = ")
      # print(nrow(X))

      # Generate a new tree based on the current
      new_trees[[j]] = update_tree(y = y_scale,
                                   X = X,
                                   type = type,
                                   curr_tree = curr_trees[[j]],
                                   node_min_size = node_min_size,
                                   s = s,
                                   max_bad_trees = max_bad_trees,
                                   splitting_rules = splitting_rules)

      # CURRENT TREE: compute the log of the marginalised likelihood + log of the tree prior
      l_old = TVPtree_full_conditional(curr_trees[[j]],
                                       Lmatleaf,
                                    current_partial_residuals,
                                    sigma2,
                                    V,
                                    inv_V,
                                    nu,
                                    lambda,
                                    tau_b,
                                    coeff_prior_conj) # + get_tree_prior(curr_trees[[j]], alpha, beta)

      # NEW TREE: compute the log of the marginalised likelihood + log of the tree prior
      l_new = TVPtree_full_conditional(new_trees[[j]],
                                       Lmatleaf,
                                    current_partial_residuals,
                                    sigma2,
                                    V,
                                    inv_V,
                                    nu,
                                    lambda,
                                    tau_b,
                                    coeff_prior_conj) # + get_tree_prior(new_trees[[j]], alpha, beta)


      # print("l_old = ")
      # print(l_old)
      #
      # print("l_new = ")
      # print(l_new)
      #
      # print("curr_trees[[j]] = ")
      # print(curr_trees[[j]])
      #
      #
      # print("new_trees[[j]] = ")
      # print(new_trees[[j]])

      # # Exponentiate the results above
      # if(type == 'grow'){
      #   # a = exp(l_new - l_old)*ratio_grow(new_trees[[j]], curr_trees[[j]])
      #   a = exp(l_new - l_old)*ratio_grow(curr_trees[[j]], new_trees[[j]])
      # } else if(type == 'prune'){
      #   # a = exp(l_new - l_old)*ratio_prune(new_trees[[j]], curr_trees[[j]])
      #   a = exp(l_new - l_old)*ratio_prune(curr_trees[[j]], new_trees[[j]])
      # } else{
      #   a = exp(l_new - l_old)
      # }


      a <- get_MH_probability2(curr_trees[[j]], new_trees[[j]], l_old, l_new,
                               type, trans_prob,
                               alpha, beta)




      if(a > runif(1)) {

        curr_trees[[j]] = new_trees[[j]] # The current tree "becomes" the new tree, if the latter is better

        if (type =='change'){
          var_count[curr_trees[[j]]$var[1] - 1] = var_count[curr_trees[[j]]$var[1] - 1] - 1
          var_count[curr_trees[[j]]$var[2] - 1] = var_count[curr_trees[[j]]$var[2] - 1] + 1
        }

        if (type=='grow'){
          var_count[curr_trees[[j]]$var - 1] = var_count[curr_trees[[j]]$var - 1] + 1 } # -1 because of the intercept in X

        if (type=='prune'){
          var_count[curr_trees[[j]]$var - 1] = var_count[curr_trees[[j]]$var - 1] - 1 } # -1 because of the intercept in X
      }

      # Update mu whether tree accepted or not
      curr_trees[[j]] = TVPsimulate_beta(curr_trees[[j]],
                                         Lmatleaf,
                                         current_partial_residuals,
                                         sigma2,
                                         inv_V,
                                         tau_b,
                                         nu,
                                         coeff_prior_conj)

      current_fit = TVPget_predictions(curr_trees[[j]], Lmatleaf, X, single_tree = TRUE)
      predictions = predictions - tree_fits_store[,j] # subtract the old fit
      predictions = predictions + current_fit # add the new fit
      tree_fits_store[,j] = current_fit # update the new fit

    } # End loop through trees

    # Updating the predictions (y_hat)
    # predictions = get_predictions(curr_trees, X, single_tree = ntrees == 1, ancestors)

    sum_of_squares = sum((y_scale - predictions)^2)

    # Update sigma2 (variance of the residuals)
    sigma2 = update_sigma2(sum_of_squares, n = length(y_scale), nu, lambda)

    # Update sigma2_beta0 and sigma2_beta1
    if (vars_inter_slope == TRUE) {
      vars_betas = TVPupdate_vars_intercepts_slopes(curr_trees, ntrees, sigma2, coeff_prior_conj)
      V = 1/c(#vars_betas$var_inter,
              vars_betas$var_slopes)
      inv_V = 1/V
    }

    # # Update s = (s_1, ..., s_p), where s_p is the probability that predictor p is used to create new terminal nodes
    # if (sparse == 'TRUE' & i > floor(TotIter*0.1)){
    #   s = update_s(var_count, p, 1)
    # }


    # Update s = (s_1, ..., s_p), where s_p is the probability that predictor q in 1:p is used to create new terminal nodes
    if ( (sparse ==TRUE) & (i > floor(nburn * 0.5))) {
      s_update <- update_s(var_count, p, alpha_s)
      s <- s_update[[1]]
      if(length(s) != ncol(X_orig)){
        print("s = ")
        print(s)
        print("length(s) = ")
        print(length(s))
        print("ncol(X) = ")
        print(ncol(X_orig))
        print("p = ")
        print(p)

      }

      if(alpha_prior == TRUE){
        alpha_s <- update_alpha(s, alpha_scale, alpha_a, alpha_b, p, s_update[[2]])
      }
    }

    # if(!is.numeric(V[1])){
    #   print("V = ")
    #   print(V)
    #   stop("!is,numeric(V[1])")
    # }


    # If at the right place, store everything
    if((i > nburn) & ((i - nburn) %% nthin) == 0) {
      curr = (i - nburn)/nthin
      tree_store[[curr]] = curr_trees
      sigma2_store[curr] = sigma2
      y_hat_store[curr,] = predictions
      var_count_store[curr,] = var_count
      s_prob_store[curr,] = s
      vars_betas_store[curr,] = V
      alpha_s_store[curr] <- alpha_s

    }


  } # End iterations loop

  cat('\n') # Make sure progress bar ends on a new line

  return(list(trees = tree_store,
              sigma2 = sigma2_store*y_sd^2,
              y_hat =(y_hat_store + (y_max + y_min)/2 )*y_sd + y_mean, # y_hat_store*y_sd + y_mean,
              center_x = center,
              scale_x = scale,
              npost = npost,
              nburn = nburn,
              nthin = nthin,
              ntrees = ntrees,
              y_mean = y_mean,
              y_sd = y_sd,
              ancestors = ancestors,
              var_count_store = var_count_store,
              s = s_prob_store,
              vars_betas = vars_betas_store,
              ntrain = n,
              y_max = y_max,
              y_min = y_min))

} # End main function






#///////////////////////////////////////////////////////////////////////////////////////////////////////////
# MOTR-BART for classification
#///////////////////////////////////////////////////////////////////////////////////////////////////////////

#' @export

motr_bart_class = function(x,
                           y,
                           sparse = TRUE,
                           vars_inter_slope = TRUE,
                           ntrees = 10,
                           node_min_size = 5,
                           alpha = 0.95,
                           beta = 2,
                           nu = 3,
                           lambda = 0.1,
                           sigma2 = 1,
                           nburn = 1000,
                           npost = 1000,
                           nthin = 1,
                           ancestors = FALSE) {

  y = as.integer(as.factor(y)) -1
  x = as.data.frame(x)

  # Quantities needed for prediction
  center = apply(x, 2, mean)
  scale = apply(x, 2, sd)
  aux.X = lapply(x, unique) # Checking how many unique values each variable has
  unique.values.X = unlist(lapply(aux.X, length))

  center[which(unique.values.X<=2)] = 0
  scale[which(unique.values.X<=2)] = 1

  X_orig = x
  X = as.matrix(cbind(1,scale(x, center, scale))) # standardising the covariates and adding an intercept

  # Extract control parameters
  node_min_size = node_min_size

  # Extract MCMC details
  TotIter = nburn + npost*nthin # Total of iterations

  # Storage containers
  store_size = npost
  tree_store = vector('list', store_size)
  sigma2_store = rep(NA, store_size)
  y_hat_store = matrix(NA, ncol = length(y), nrow = store_size)
  var_count = rep(0, ncol(X_orig))
  var_count_store = matrix(0, ncol = ncol(X_orig), nrow = store_size)
  s_prob_store = matrix(0, ncol = ncol(X_orig), nrow = store_size)
  tree_fits_store = matrix(0, ncol = ntrees, nrow = length(y))

  # Scale the response target variable
  y_mean = mean(y)
  y_sd = sd(y)
  y_scale = (y - y_mean)/y_sd
  n = length(y_scale)
  p = ncol(X_orig)
  s = rep(1/p, p)

  # Prior for the beta vector
  tau_b = ntrees
  V = rep(1/tau_b, 2)
  inv_V = 1/V

  # Initial values
  z = ifelse(y == 0, -3, 3)

  # Create a list of trees for the initial stump
  curr_trees = create_stump(num_trees = ntrees,
                            y = z,
                            X = X)
  # Initialise the new trees as current one
  new_trees = curr_trees

  # Initialise the predicted values to zero
  predictions = get_predictions(curr_trees, X, single_tree = ntrees == 1, ancestors)

  # Set up a progress bar
  pb = utils::txtProgressBar(min = 1, max = TotIter,
                             style = 3, width = 60,
                             title = 'Running rBART...')

  # Start the MCMC iterations loop
  for (i in 1:TotIter) {

    utils::setTxtProgressBar(pb, i)

    # If at the right place, store everything
    if((i > nburn) & ((i - nburn) %% nthin) == 0) {
      curr = (i - nburn)/nthin
      tree_store[[curr]] = curr_trees
      sigma2_store[curr] = sigma2
      y_hat_store[curr,] = pnorm(predictions)
      var_count_store[curr,] = var_count
      s_prob_store[curr,] = s
    }

    # Start looping through trees
    for (j in 1:ntrees) {

      # Calculate partial residuals for current tree
      # if(ntrees > 1) {
      #   current_partial_residuals = z -
      #     get_predictions(curr_trees[-j], X, single_tree = ntrees == 2, ancestors)
      # } else {
      #   current_partial_residuals = z
      # }

      current_partial_residuals = z - predictions + tree_fits_store[,j]

      # Propose a new tree via grow/change/prune/swap
      type = sample(c('grow', 'prune', 'change', 'swap'), 1)
      if(i < max(floor(0.1*nburn), 10)) type = 'grow' # Grow for the first few iterations

      # Generate a new tree based on the current
      new_trees[[j]] = update_tree(y = z,
                                   X = X,
                                   type = type,
                                   curr_tree = curr_trees[[j]],
                                   node_min_size = node_min_size,
                                   s)

      # NEW TREE: compute the log of the marginalised likelihood + log of the tree prior
      l_new = tree_full_conditional(new_trees[[j]],
                                    X,
                                    current_partial_residuals,
                                    sigma2,
                                    V,
                                    inv_V,
                                    nu,
                                    lambda,
                                    tau_b,
                                    ancestors) +
        get_tree_prior(new_trees[[j]], alpha, beta)

      # CURRENT TREE: compute the log of the marginalised likelihood + log of the tree prior
      l_old = tree_full_conditional(curr_trees[[j]],
                                    X,
                                    current_partial_residuals,
                                    sigma2,
                                    V,
                                    inv_V,
                                    nu,
                                    lambda,
                                    tau_b,
                                    ancestors) +
        get_tree_prior(curr_trees[[j]], alpha, beta)

      # Exponentiate the results above
      a = exp(l_new - l_old)

      # The current tree "becomes" the new tree, if the latter is better
      if(a > runif(1)) {
        curr_trees[[j]] = new_trees[[j]]

        if (type =='change'){
          var_count[curr_trees[[j]]$var[1] - 1] = var_count[curr_trees[[j]]$var[1] - 1] - 1
          var_count[curr_trees[[j]]$var[2] - 1] = var_count[curr_trees[[j]]$var[2] - 1] + 1
        }

        if (type=='grow'){
          var_count[curr_trees[[j]]$var - 1] = var_count[curr_trees[[j]]$var - 1] + 1 } # -1 because of the intercept in X

        if (type=='prune'){
          var_count[curr_trees[[j]]$var - 1] = var_count[curr_trees[[j]]$var - 1] - 1 } # -1 because of the intercept in X
      }

      # Update mu whether tree accepted or not
      curr_trees[[j]] = simulate_beta(curr_trees[[j]],
                                      X,
                                      current_partial_residuals,
                                      sigma2,
                                      inv_V,
                                      tau_b,
                                      nu,
                                      ancestors)

      current_fit = get_predictions(curr_trees[j], X, single_tree = TRUE, ancestors)
      predictions = predictions - tree_fits_store[,j] # subtract the old fit
      predictions = predictions + current_fit # add the new fit
      tree_fits_store[,j] = current_fit # update the new fit

    } # End loop through trees

    # Updating the predictions (y_hat)
    # predictions = get_predictions(curr_trees, X, single_tree = ntrees == 1, ancestors)

    # Update z (latent variable)
    z = update_z(y, predictions)

    # sum_of_squares = sum((y_scale - pnorm(predictions))^2)

    # Update sigma2_beta0 and sigma2_beta1
    if (vars_inter_slope == 'TRUE') {
      vars_betas = update_vars_intercepts_slopes(curr_trees, ntrees, sigma2)
      V = 1/c(vars_betas$var_inter, vars_betas$var_slopes)
      inv_V = 1/V
    }

    # Update sigma2 (variance of the residuals)
    # sigma2 = update_sigma2(sum_of_squares, n = n, nu, lambda)

    # Update s = (s_1, ..., s_p), where s_p is the probability that predictor p is used to create new terminal nodes
    if (sparse == 'TRUE' & i > floor(TotIter*0.1)){
      s = update_s(var_count, p, 1)
    }

  } # End iterations loop

  cat('\n') # Make sure progress bar ends on a new line

  return(list(trees = tree_store,
              sigma2 = sigma2_store,
              y_hat = y_hat_store,
              center_x = center,
              scale_x = scale,
              npost = npost,
              nburn = nburn,
              nthin = nthin,
              ntrees = ntrees,
              y_mean = y_mean,
              y_sd = y_sd,
              ancestors = ancestors,
              var_count_store = var_count_store,
              s = s_prob_store))

} # End main function
