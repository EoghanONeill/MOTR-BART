# -------------------------------------------------------------------------#
# Description: this script contains auxiliar functions needed to update    #
# the trees with details and to map the predicted values to each obs       #
# -------------------------------------------------------------------------#

# 1. fill_tree_details: takes a tree matrix and returns the number of obs in each node in it and the indices of each observation in each terminal node
# 2. get_predictions: gets the predicted values from a current set of trees
# 3. get_children: it's a function that takes a node and, if the node is terminal, returns the node. If not, returns the children and calls the function again on the children
# 4. resample: an auxiliar function
# 5. get_ancestors: get the ancestors of all terminal nodes in a tree
# 6. update_s: full conditional of the vector of splitting probability.
# 7. update_vars_intercepts_slopes: updates the variances of the intercepts and slopes

# Fill_tree_details -------------------------------------------------------

fill_tree_details = function(curr_tree, X) {

  # Collect right bits of tree
  tree_matrix = curr_tree$tree_matrix

  # Create a new tree matrix to overwrite
  new_tree_matrix = tree_matrix

  # Start with dummy node indices
  node_indices = rep(1, nrow(X))

  if (nrow(tree_matrix) > 1) {

    # For all but the top row, find the number of observations falling into each one
    for(i in 2:nrow(tree_matrix)) {

      # Get the parent
      curr_parent = as.numeric(tree_matrix[i,'parent'])

      # Find the split variable and value of the parent
      split_var = as.numeric(tree_matrix[curr_parent,'split_variable'])
      split_val = as.numeric(tree_matrix[curr_parent, 'split_value'])

      # # Find whether it's a left or right terminal node
      # left_or_right = ifelse(tree_matrix[curr_parent,'child_left'] == i,
      #                        'left', 'right')
      # if(left_or_right == 'left') {
      if (tree_matrix[curr_parent, "child_left"] == i) {
        # If left use less than condition
        new_tree_matrix[i,'node_size'] <- sum(X[node_indices == curr_parent,split_var] < split_val)
        node_indices[node_indices == curr_parent][X[node_indices == curr_parent,split_var] < split_val] <- i
      } else {
        # # If right use greater than condition
        # new_tree_matrix[i,'node_size'] = sum(X[node_indices == curr_parent,split_var] >= split_val)
        # node_indices[node_indices == curr_parent][X[node_indices == curr_parent,split_var] >= split_val] = i

        new_tree_matrix[i, "node_size"] <- sum(node_indices == curr_parent)
        node_indices[node_indices == curr_parent] <- i
      }
    } # End of loop through table
  }


  return(list(tree_matrix = new_tree_matrix,
              node_indices = node_indices))

} # End of function

# Get predictions ---------------------------------------------------------

get_predictions = function(trees, X, single_tree = FALSE, ancestors) {

  # Stop nesting problems in case of multiple trees
  if(is.null(names(trees)) & (length(trees) == 1)) trees = trees[[1]]

  # Normally trees will be a list of lists but just in case
  if(single_tree) {

    # Deal with just a single tree
    if(nrow(trees$tree_matrix) == 1) {
      # predictions = rep(trees$tree_matrix[1, 'mu'], nrow(X))
      beta_hat = as.numeric(unlist(strsplit(trees$tree_matrix[1, 'beta_hat'],",")))
      predictions = rep(beta_hat[1], nrow(X))

    } else {
      # Loop through the node indices to get predictions
      predictions = rep(NA, nrow(X))
      unique_node_indices = unique(trees$node_indices)
      # Get the node indices for the current X matrix
      curr_X_node_indices = fill_tree_details(trees, X)$node_indices

      # if (ancestors == FALSE) {lm_vars <- c(1, sort(unique(as.numeric(split_vars_tree))))}
      if (ancestors == FALSE) {
        which_internal = which(trees$tree_matrix[,'terminal'] == 0)
        split_vars_tree <- trees$tree_matrix[which_internal, 'split_variable']
        lm_vars <- c(1, sort_unique(as.numeric(split_vars_tree)))
      }
      #if (ancestors == 'all covariates') {lm_vars <- 1:ncol(X)}
      if (ancestors == TRUE) {get_ancs <- get_ancestors(trees)}

      n = nrow(X)

      # Now loop through all node indices to fill in details
      for(i in 1:length(unique_node_indices)) {
        if (ancestors == TRUE) {
          lm_vars = c(1, get_ancs[which(get_ancs[,'terminal'] == unique_node_indices[i]), 'ancestor']) # Get the corresponding ancestors of the current terminal node
        }
        X_node = matrix(X[,lm_vars], nrow=n)[curr_X_node_indices == unique_node_indices[i],]
        beta_hat = as.numeric(unlist(strsplit(trees$tree_matrix[unique_node_indices[i], 'beta_hat'],",")))
        predictions[curr_X_node_indices == unique_node_indices[i]] = X_node%*%beta_hat
      }
    }
    # More here to deal with more complicated trees - i.e. multiple trees
  } else {

    # Do a recursive call to the function
    partial_trees = trees
    partial_trees[[1]] = NULL # Blank out that element of the list
    predictions = get_predictions(trees[[1]], X, single_tree = TRUE, ancestors)  +
      get_predictions(partial_trees, X,
                      single_tree = length(partial_trees) == 1, ancestors)
    #single_tree = !is.null(names(partial_trees)))
    # The above only sets single_tree to if the names of the object is not null (i.e. is a list of lists)
  }

  return(predictions)
}





TVPget_predictions = function(trees, Lmat, X,  single_tree = FALSE) {

  # Stop nesting problems in case of multiple trees
  if(is.null(names(trees)) & (length(trees) == 1)) trees = trees[[1]]

  # Normally trees will be a list of lists but just in case
  if(single_tree) {

    # Deal with just a single tree
    if(nrow(trees$tree_matrix) == 1) {
      # predictions = rep(trees$tree_matrix[1, 'mu'], nrow(X))
      beta_hat = as.numeric(unlist(strsplit(trees$tree_matrix[1, 'beta_hat'],",")))
      # predictions = rep(beta_hat[1], nrow(X))
      predictions <- Lmat%*%beta_hat

    } else {
      # Loop through the node indices to get predictions
      predictions = rep(NA, nrow(X))
      unique_node_indices = unique(trees$node_indices)
      # Get the node indices for the current X matrix
      curr_X_node_indices = fill_tree_details(trees, X)$node_indices
      # which_internal = which(trees$tree_matrix[,'terminal'] == 0)
      # split_vars_tree <- trees$tree_matrix[which_internal, 'split_variable']

      # # if (ancestors == FALSE) {lm_vars <- c(1, sort(unique(as.numeric(split_vars_tree))))}
      # if (ancestors == FALSE) {lm_vars <- c(1, sort_unique(as.numeric(split_vars_tree)))}
      # #if (ancestors == 'all covariates') {lm_vars <- 1:ncol(X)}
      # if (ancestors == TRUE) {get_ancs <- get_ancestors(trees)}

      # n = nrow(X)

      # Now loop through all node indices to fill in details
      for(i in 1:length(unique_node_indices)) {
        # if (ancestors == TRUE) {
        #   lm_vars = c(1, get_ancs[which(get_ancs[,'terminal'] == unique_node_indices[i]), 'ancestor']) # Get the corresponding ancestors of the current terminal node
        # }
        # X_node = Lmat[curr_X_node_indices == unique_node_indices[i],]
        beta_hat = as.numeric(unlist(strsplit(trees$tree_matrix[unique_node_indices[i], 'beta_hat'],",")))
        predictions[curr_X_node_indices == unique_node_indices[i]] =  Lmat[curr_X_node_indices == unique_node_indices[i],]%*%beta_hat
      }
    }
    # More here to deal with more complicated trees - i.e. multiple trees
  } else {

    # Do a recursive call to the function
    partial_trees = trees
    partial_trees[[1]] = NULL # Blank out that element of the list
    predictions = TVPget_predictions(trees[[1]], Lmat, X, single_tree = TRUE)  +
      TVPget_predictions(partial_trees, Lmat, X,
                      single_tree = length(partial_trees) == 1)
    #single_tree = !is.null(names(partial_trees)))
    # The above only sets single_tree to if the names of the object is not null (i.e. is a list of lists)
  }

  return(predictions)
}

# get_children ------------------------------------------------------------

get_children = function(tree_mat, parent) {
  # Create a holder for the children
  all_children = NULL
  if(as.numeric(tree_mat[parent,'terminal']) == 1) {
    # If the node is terminal return the list so far
    return(c(all_children, parent))
  } else {
    # If not get the current children
    curr_child_left = as.numeric(tree_mat[parent, 'child_left'])
    curr_child_right = as.numeric(tree_mat[parent, 'child_right'])
    # Return the children and also the children of the children recursively
    return(c(all_children,
             get_children(tree_mat,curr_child_left),
             get_children(tree_mat,curr_child_right)))
  }
}

# Sample function ----------------------------------------------------------

resample <- function(x, ...) x[sample.int(length(x), size=1), ...]

# Get ancestors of each terminal node --------------------------------------

get_ancestors = function(tree){

  save_ancestor = NULL
  which_terminal = which(tree$tree_matrix[,'terminal'] == 1)

  if(nrow(tree$tree_matrix) == 1) {
    save_ancestor = cbind(terminal = NULL,
                          ancestor = NULL)
  } else {
    for (k in 1:length(which_terminal)){
      get_parent = as.numeric(as.character(tree$tree_matrix[which_terminal[k], 'parent'])) # get the 1st parent
      get_split_var = as.numeric(as.character(tree$tree_matrix[get_parent, 'split_variable'])) # then, get the covariate associated to the row of the parent

      save_ancestor = rbind(save_ancestor,
                            cbind(terminal = which_terminal[k],
                                  # parent   = get_parent,
                                  ancestor = get_split_var))
      while (get_parent > 1){
        get_parent = as.numeric(as.character(tree$tree_matrix[get_parent,'parent'])) # then, get the subsequent parent
        get_split_var = as.numeric(as.character(tree$tree_matrix[get_parent, 'split_variable'])) # then, get the covariate associated to the row of the new parent
        save_ancestor = rbind(save_ancestor,
                              cbind(terminal = which_terminal[k],
                                    # parent   = get_parent,
                                    ancestor = get_split_var))
      }
    }
    save_ancestor = unique(save_ancestor) # remove duplicates
    save_ancestor = save_ancestor[order(save_ancestor[,1], save_ancestor[,2]),] # sort by terminal and ancestor
  }

  return(save_ancestor)
}

# update_s = function(var_count, p, alpha_s){
#   s_ = rdirichlet(1, alpha_s/p + var_count)
#   return(s_)
# }


update_s <- function(var_count, p, alpha_s) {
  # s_ = rdirichlet(1, as.vector((alpha_s / p ) + var_count))

  # // Get shape vector
  # shape_up = alpha_s / p
  shape_up = as.vector((alpha_s / p ) + var_count)

  # // Sample unnormalized s on the log scale
  templogs = rep(NA, p)
  for(i in 1:p) {
    templogs[i] = SoftBart:::rlgam(shape = shape_up[i])
  }

  if(any(templogs== -Inf)){
    print("alpha_s = ")
    print(alpha_s)
    print("var_count = ")
    print(var_count)
    print("templogs = ")
    print(templogs)
    stop('templogs == -Inf')
  }

  # // Normalize s on the log scale, then exponentiate
  # templogs = templogs - log_sum_exp(hypers.logs);
  max_log = max(templogs)
  templogs2 = templogs - (max_log + log(sum(exp( templogs  -  max_log ))))


  s_ = exp(templogs2)

  # if(any(s_==0)){
  #   print("templogs2 = ")
  #   print(templogs2)
  #   print("templogs = ")
  #   print(templogs)
  #   print("alpha_s = ")
  #   print(alpha_s)
  #   print("var_count = ")
  #   print(var_count)
  #   print("s_ = ")
  #   print(s_)
  #   stop('s_ == 0')
  # }

  ret_list <- list()
  ret_list[[1]] <- s_
  ret_list[[2]] <- mean(templogs2)


  return(ret_list)
}





update_alpha <- function(s, alpha_scale, alpha_a, alpha_b, p, mean_log_s) {

  # create inputs for likelihood

  # log_s <- log(s)
  # mean_log_s <- mean(log_s)
  # p <- length(s)
  # alpha_scale   # denoted by lambda_a in JRSSB paper

  rho_grid <- (1:1000)/1001

  alpha_grid <- alpha_scale * rho_grid / (1 - rho_grid )

  logliks <- alpha_grid * mean_log_s +
    lgamma(alpha_grid) -
    p*lgamma(alpha_grid/p) +
    (alpha_a - 1)*log(rho_grid) + (alpha_b-1)*log(1- rho_grid)
  # dbeta(x = rho_grid, shape1 = alpha_a, shape2 = alpha_b, ncp = 0, log = TRUE)


  # logliks <- log(ddirichlet( t(matrix(s, p, 1000))  , t(matrix( rep(alpha_grid/p,p) , p , 1000)  ) ) ) +
  #   (alpha_a - 1)*log(rho_grid) + (alpha_b-1)*log(1- rho_grid)
  # # dbeta(x = rho_grid, shape1 = alpha_a, shape2 = alpha_b, ncp = 0, log = TRUE)

  # logliks <- rep(NA, 1000)
  # for(i in 1:1000){
  #   logliks[i] <- log(ddirichlet(s  , rep(alpha_grid[i]/p,p) ) ) +
  #     (alpha_a - 1)*log(rho_grid[i]) + (alpha_b-1)*log(1- rho_grid[i])
  # }

  max_ll <- max(logliks)
  logsumexps <- max_ll + log(sum(exp( logliks  -  max_ll )))

  # print("logsumexps = ")
  # print(logsumexps)

  logliks <- exp(logliks - logsumexps)

  if(any(is.na(logliks))){
    print("logliks = ")
    print(logliks)

    print("logsumexps = ")
    print(logsumexps)

    print("mean_log_s = ")
    print(mean_log_s)

    print("lgamma(alpha_grid) = ")
    print(lgamma(alpha_grid))

    print("p*lgamma(alpha_grid/p) = ")
    print(p*lgamma(alpha_grid/p))

    print("(alpha_a - 1)*log(rho_grid) + (alpha_b-1)*log(1- rho_grid) = ")
    print((alpha_a - 1)*log(rho_grid) + (alpha_b-1)*log(1- rho_grid))

    print("max_ll = ")
    print(max_ll)

    # print("s = ")
    # print(s)


  }

  # print("logliks = ")
  # print(logliks)

  rho_ind <- sample.int(1000,size = 1, prob = logliks)


  return(alpha_grid[rho_ind])
}





update_vars_intercepts_slopes <- function(trees, n_tress, sigma2, coeff_prior_conj, a0 = 1, b0 = 1, a1 = 1, b1 = 1){

  n_terminal = 0
  n_vars_terminal = 0
  sum_of_squares_inter = 0
  sum_of_squares_slopes = 0

  for (i in 1:n_tress) {
    # Get current set of trees
    tree = trees[[i]]
    # get the terminal nodes
    terminal_nodes = as.numeric(which(tree$tree_matrix[,'terminal'] == 1))
    # get all coefficients of the linear predictors for each terminal node
    all_coef = strsplit(tree$tree_matrix[terminal_nodes, 'beta_hat'], ',')
    # get intercepts
    inter = as.numeric(unlist(lapply(all_coef, '[', 1)))
    # get slopes
    slopes = as.numeric(unlist(lapply(all_coef, '[', -1)))

    n_terminal = n_terminal + length(terminal_nodes)
    n_vars_terminal = n_vars_terminal + length(slopes)
    sum_of_squares_inter = sum_of_squares_inter + sum(inter^2)
    sum_of_squares_slopes = sum_of_squares_slopes + sum(slopes^2)
  }

  if(coeff_prior_conj == TRUE){
    return(list(var_inter = rgamma(1, (n_terminal/2) + a0, sum_of_squares_inter/(2*sigma2) + b0),
                var_slopes = rgamma(1, (n_vars_terminal/2) + a1, sum_of_squares_slopes/(2*sigma2) + b1)))
  }else{
    return(list(var_inter = rgamma(1, (n_terminal/2) + a0, sum_of_squares_inter/(2) + b0),
                var_slopes = rgamma(1, (n_vars_terminal/2) + a1, sum_of_squares_slopes/(2) + b1)))
  }


}



TVPupdate_vars_intercepts_slopes <- function(trees, n_tress, sigma2, coeff_prior_conj, a0 = 1, b0 = 1, a1 = 1, b1 = 1){

  n_terminal = 0
  n_vars_terminal = 0
  # sum_of_squares_inter = 0
  sum_of_squares_slopes = 0

  for (i in 1:n_tress) {
    # Get current set of trees
    tree = trees[[i]]
    # get the terminal nodes
    terminal_nodes = as.numeric(which(tree$tree_matrix[,'terminal'] == 1))
    # get all coefficients of the linear predictors for each terminal node
    all_coef = strsplit(tree$tree_matrix[terminal_nodes, 'beta_hat'], ',')
    # get intercepts
    # inter = as.numeric(unlist(lapply(all_coef, '[', 1)))
    # get slopes
    # slopes = as.numeric(unlist(lapply(all_coef, '[', -1)))
    slopes = as.numeric(unlist(all_coef))

    # n_terminal = n_terminal + length(terminal_nodes)
    n_vars_terminal = n_vars_terminal + length(slopes)
    # sum_of_squares_inter = sum_of_squares_inter + sum(inter^2)
    sum_of_squares_slopes = sum_of_squares_slopes + sum(slopes^2)
  }
  if(coeff_prior_conj == TRUE){
    return(list(#var_inter = rgamma(1, (n_terminal/2) + a0, sum_of_squares_inter/(2*sigma2) + b0),
                var_slopes = rgamma(1, (n_vars_terminal/2) + a1, sum_of_squares_slopes/(2*sigma2) + b1)))
  }else{
    return(list(#var_inter = rgamma(1, (n_terminal/2) + a0, sum_of_squares_inter/(2*sigma2) + b0),
      var_slopes = rgamma(1, (n_vars_terminal/2) + a1, sum_of_squares_slopes/(2) + b1)))
  }
}




sample_move = function(curr_tree, i, nburn, trans_prob){

  if (nrow(curr_tree$tree_matrix) == 1 || i < max(floor(0.1*nburn), 5)) {
    type = 'grow'
  } else {
    type = sample(c('grow', 'prune', 'change'),  1, prob = trans_prob)
  }
  return(type)
}


# This function returns the Metropolis-Hastings acceptance probability in line
# with Kapelner, A., and Bleich, J. (2013). "bartMachine: Machine learning with
# Bayesian additive regression trees."
get_MH_probability2 <- function(curr_tree, new_tree,l_old, l_new,
                                type, trans_prob,
                                alpha, beta) {
  # Number of terminal nodes in current tree
  # b_j <- sum(curr_tree$tree_matrix[, "terminal"])

  # Get the tree type probabilities
  # if(nrow(curr_tree$tree_matrix) == 1 ){
  #   prob_grow <-  1 #trans_prob[1]
  #   prob_prune <- trans_prob[2]
  # }else{
  #   prob_grow <- trans_prob[1]
  #   prob_prune <- trans_prob[2]
  # }

  prob_grow <- trans_prob[1]
  prob_prune <- trans_prob[2]
  # l_new <- get_logL(new_tree, new_partial_resid_rescaled, att_weights_new, mu_mu, sigma2_mu, sigma2)
  # l_old <- get_logL(curr_tree, curr_partial_resid_rescaled, att_weights_current, mu_mu, sigma2_mu, sigma2)

  if(type == 'grow'){
    # a = exp(l_new - l_old + get_tree_prior(new_tree, alpha, beta)  - get_tree_prior(curr_tree, alpha, beta) )*
    #   ratio_grow(curr_tree, new_tree) * (prob_prune / prob_grow)
    a = exp(l_new - l_old + get_tree_prior(new_tree, alpha, beta)  - get_tree_prior(curr_tree, alpha, beta) +
              log(ratio_grow(curr_tree, new_tree ) ))*
      (prob_prune / prob_grow)
  } else if(type == 'prune'){
    a = exp(l_new - l_old + get_tree_prior(new_tree, alpha, beta)  - get_tree_prior(curr_tree, alpha, beta) +
              log( ratio_prune(curr_tree, new_tree) ))*
      (prob_grow / prob_prune)
  } else{
    a = exp(l_new - l_old)
  }

  if(is.na(a)){
    print("get_tree_prior(new_tree, alpha, beta) = ")
    print(get_tree_prior(new_tree, alpha, beta))

    print("get_tree_prior(curr_tree, alpha, beta) ) = ")
    print(get_tree_prior(curr_tree, alpha, beta) )


    print("ratio_grow(curr_tree, new_tree ) ) = ")
    print(ratio_grow(curr_tree, new_tree ) )

    print("ratio_prune(curr_tree, new_tree)  = ")
    print(ratio_prune(curr_tree, new_tree) )

    print("new_tree = ")
    print(new_tree)

    print("curr_tree = ")
    print(curr_tree)

    print("alpha = ")
    print(alpha)

    print("beta = ")
    print(beta)

    print("l_new = ")
    print(l_new)

    print("l_old = ")
    print(l_old)


  }

  # r <- exp(l_new - l_old) * transition_ratio * tree_ratio
  return(min(1, a))
}


# This function returns the number of second generation nodes "w"
get_gen2 <- function(tree) {
  if (nrow(tree$tree_matrix) == 1) {
    w <- 0
  } else {
    # indeces <- which(tree$tree_matrix[, "terminal"] == 1)
    # Determine the parent for each terminal node and find the duplicated parents
    # w <- as.numeric(sum(duplicated(tree$tree_matrix[indeces,'parent'])))
    # parents <- tree$tree_matrix[indeces, "parent"]
    parents <- tree$tree_matrix[tree$tree_matrix[, "terminal"] == 1, "parent"]
    w <- parents[duplicated(parents)]
  }
  return(w)
}
