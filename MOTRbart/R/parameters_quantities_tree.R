# -------------------------------------------------------------------------#
# Description: this script contains 2 functions that are used to generate  #
#              the predictions, update variance and compute the tree prior #
#              and the marginalised likelihood                             #
# -------------------------------------------------------------------------#

# 1. simulate_mu: generate the 'betas' in the linear predictors
# 2. updata_sigma2: updates the parameters sigma2
# 3. get_tree_prior: returns the tree log prior score
# 4. tree_full_conditional: computes the marginalised likelihood for all nodes for a given tree

# Compute the full conditionals -------------------------------------------------

tree_full_conditional = function(tree, X, R, sigma2, V, inv_V, nu, lambda, tau_b) {

  # Select the lines that correspond to terminal and internal nodes
  which_terminal = which(tree$tree_matrix[,'terminal'] == 1)
  which_internal = which(tree$tree_matrix[,'terminal'] == 0)

  # Get the node indices for each terminal node
  curr_X_node_indices = tree$node_indices
  unique_node_indices = unique(tree$node_indices)
  log_post = NULL

  # Get the covariates that have been used as a split
  split_vars_tree <- tree$tree_matrix[which_internal, 'split_variable']
  lm_vars <- c(1, sort(unique(as.numeric(split_vars_tree))))
  p = length(lm_vars)
  inv_V = diag(p)*inv_V

  # Compute the log marginalised likelihood for each terminal node
  for(i in 1:length(unique_node_indices)) {
    X_node = X[curr_X_node_indices == unique_node_indices[i], lm_vars]
    r_node = R[curr_X_node_indices == unique_node_indices[i]]
    Lambda_node_inv = t(X_node)%*%X_node + inv_V
    Lambda_node = solve(t(X_node)%*%X_node + inv_V)
    mu_node = Lambda_node%*%((t(X_node))%*%r_node)

    log_post[i] = -0.5 * log(V) +
      0.5*log(1/det(Lambda_node_inv)) -
      (1/(2*sigma2)) * (- t(mu_node)%*%Lambda_node_inv%*%mu_node)

  }
  return(sum(log_post))
}


# Simulate_par -------------------------------------------------------------

simulate_mu = function(tree, X, R, sigma2, inv_V, tau_b, nu) {

  # First find which rows are terminal and internal nodes
  which_terminal = which(tree$tree_matrix[,'terminal'] == 1)
  which_internal = which(tree$tree_matrix[,'terminal'] == 0)

  # Get node indices
  curr_X_node_indices = tree$node_indices
  unique_node_indices = unique(tree$node_indices)

  # Wipe all the old parameters out for other nodes
  tree$tree_matrix[,'mu'] = NA

  # Get the covariates that have been used as a split
  split_vars_tree <- tree$tree_matrix[which_internal, 'split_variable']
  lm_vars <- c(1, sort(unique((as.numeric(split_vars_tree)))))
  p = length(lm_vars)
  inv_V = diag(p)*inv_V

  for(i in 1:length(unique_node_indices)) {
    X_node = X[curr_X_node_indices == unique_node_indices[i], lm_vars] # Only variables that have been used as split
    r_node = R[curr_X_node_indices == unique_node_indices[i]]
    Lambda_node = solve(t(X_node)%*%X_node + inv_V)

    # Generate betas  -------------------------------------------------
    mu = rmvnorm(1,
                 mean = Lambda_node%*%(t(X_node)%*%r_node),
                 sigma = sigma2*Lambda_node)

    # Put in just the ones that are useful, otherwise 0.
    aux_mu = rep(0, ncol(X))
    aux_mu[lm_vars] = mu # Only variables that have been used as split
    tree$tree_matrix[unique_node_indices[i],'mu'] = paste(aux_mu, collapse = ',')
  }

  return(tree)
}

# Update sigma2 -------------------------------------------------------------

update_sigma2 <- function(S, n, nu, lambda){
  u = 1/rgamma(1, shape = (n + nu)/2, rate = (S + nu*lambda)/2)
  return(u)
}

# Get tree priors ---------------------------------------------------------

get_tree_prior = function(tree, alpha, beta) {

  # Need to work out the depth of the tree
  # First find the level of each node, then the depth is the maximum of the level
  level = rep(NA, nrow(tree$tree_matrix))
  level[1] = 0 # First row always level 0

  # Escpae quickly if tree is just a stump
  if(nrow(tree$tree_matrix) == 1) {
    return(log(1 - alpha)) # Tree depth is 0
  }

  for(i in 2:nrow(tree$tree_matrix)) {
    # Find the current parent
    curr_parent = as.numeric(tree$tree_matrix[i,'parent'])
    # This child must have a level one greater than it's current parent
    level[i] = level[curr_parent] + 1
  }

  # Only compute for the internal nodes
  internal_nodes = which(as.numeric(tree$tree_matrix[,'terminal']) == 0)
  log_prior = 0
  for(i in 1:length(internal_nodes)) {
    log_prior = log_prior + log(alpha) - beta * log(1 + level[internal_nodes[i]])
  }
  # Now add on terminal nodes
  terminal_nodes = which(as.numeric(tree$tree_matrix[,'terminal']) == 1)
  for(i in 1:length(terminal_nodes)) {
    log_prior = log_prior + log(1 - alpha * ((1 + level[terminal_nodes[i]])^(-beta)))
  }


  return(log_prior)

}