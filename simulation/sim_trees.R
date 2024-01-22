sim_tree_dists <- function(t) {
  if(!require("ape",character.only = TRUE, quietly = TRUE)) install.packages("ape")
  tree <- ape::rtree(t,tip.label=1:t)
  
  inv_dist <- tree_dist(tree)
  return(inv_dist)
}

sim_tree <- function(t) {
  if(!require("ape",character.only = TRUE, quietly = TRUE)) install.packages("ape")
  return(ape::rtree(t,tip.label=1:t))
}

tree_dist <- function(tree) {
  if(!require("ape",character.only = TRUE, quietly = TRUE)) install.packages("ape")
  inv_dist <- 1/ape::cophenetic.phylo(tree)
  diag(inv_dist) <- 0
  
  
  inv_dist <- inv_dist[order(as.integer(rownames(inv_dist))),order(as.integer(colnames(inv_dist)))]
  return(inv_dist)
}

sim_tree_string <- function(t)
{
  tree <- sim_tree(t)
  return(list('tree_str'=write.tree(tree),'tree'=tree))
}

