sim_tree_dists <- function(k,chosen) {
  if(!require("ape",character.only = TRUE, quietly = TRUE)) quiet(install.packages("ape"),all=TRUE)
  tree <- rtree(k,tip.label=chosen)
  inv_dist <- 1/cophenetic.phylo(tree)
  diag(inv_dist) <- 0
  
  ret_mat <- matrix(c(rep(0,400)),nrow=20)
  colnames(ret_mat) <- 1:20; rownames(ret_mat) <- 1:20
  
  inv_dist <- inv_dist[order(as.integer(rownames(inv_dist))),order(as.integer(colnames(inv_dist)))]
  ret_mat[chosen,chosen] <- inv_dist
  return(ret_mat)
}