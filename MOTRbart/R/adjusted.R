# This code follows the structure of transition probabilities that is in accordance with the bartMachine paper and soft BART
# bartMachine: Machine Learning with Bayesian Additive Regression Trees by Adam Kapelner and Justin Bleich
# soft BART: Bayesian regression tree ensembles that adapt to smoothness and sparsity by Antonio R. Linero and Yun Yang

# This functions calculates the number of terminal nodes
get_nterminal = function(tree){

  indeces = which(tree$tree_matrix[,'terminal'] == '1') # determine which indeces have a terminal node label
  b = as.numeric(length(indeces)) # take the length of these indeces, to determine the number of terminal nodes/leaves
  return(b)
}

# This function calculates the number of parents with two terminal nodes/
# second generartion internal nodes as formulated in bartMachine
get_w = function(tree){
  indeces = which(tree$tree_matrix[,'terminal'] == '1') #determine which indeces have a terminal node label
  # determine the parent for each terminal node and sum the number of duplicated parents
  w = as.numeric(sum(duplicated(tree$tree_matrix[indeces,'parent'])))
  return(w)
}

# These functions calculate the grow and prune ratios respectively according to the sof BART paper
# ratio_grow = function(new_tree, cur_tree){
ratio_grow = function(curr_tree, new_tree){
    # grow_ratio = get_nterminal(cur_tree)/(get_w(new_tree)+1)
  grow_ratio = get_nterminal(curr_tree)/(get_w(new_tree)) # (get_w(new_tree)+1)

  return(as.numeric(grow_ratio))
}

# ratio_prune = function(new_tree, cur_tree){
ratio_prune = function(curr_tree, new_tree){
  # prune_ratio = get_w(new_tree)/(get_nterminal(cur_tree)-1)
  prune_ratio = get_w(curr_tree)/(get_nterminal(curr_tree)-1) #(get_nterminal(curr_tree)-1)
  return(as.numeric(prune_ratio))
}



