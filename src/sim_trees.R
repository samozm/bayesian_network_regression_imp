sim_tree_dists <- function(t) {
  if(!require("ape",character.only = TRUE, quietly = TRUE)) install.packages("ape")
  tree <- rtree(t,tip.label=1:t)
  
  inv_dist <- tree_dist(tree)
  return(inv_dist)
}

sim_tree <- function(t) {
  if(!require("ape",character.only = TRUE, quietly = TRUE)) install.packages("ape")
  return(rtree(t,tip.label=1:t))
}

tree_dist <- function(tree) {
  inv_dist <- 1/cophenetic.phylo(tree)
  diag(inv_dist) <- 0
  
  
  inv_dist <- inv_dist[order(as.integer(rownames(inv_dist))),order(as.integer(colnames(inv_dist)))]
  return(inv_dist)
}

sim_tree_string <- function(t)
{
  tree <- sim_tree(t)
  return(write.tree(tree))
}

