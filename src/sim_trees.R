sim_tree_dists <- function(k) {
  if(!require("ape",character.only = TRUE, quietly = TRUE)) install.packages("ape")
  tree <- rtree(k,tip.label=1:k)
  inv_dist <- 1/cophenetic.phylo(tree)
  diag(inv_dist) <- 0
  
  
  inv_dist <- inv_dist[order(as.integer(rownames(inv_dist))),order(as.integer(colnames(inv_dist)))]
  return(inv_dist)
}