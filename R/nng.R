#' Compute the nearest-neighbors graph + pairwise distances.
#'
#' This is just a wrapper around C++/Rcpp function.
#'
#' @param data data.frame or matrix.
#' @param k Number of neighbors at each point.
#' @param mutual Whether or not to build mutual kNN graph.
#'
#' @return A list with $graph and $distances.
#'
#' @export
nng <- function(data, k=ceil(sqrt(nrow(data))), mutual=TRUE)
{
  mdata <- data
  if (is.data.frame(data)) mdata <- as.matrix(data)
  res <- findNeighbors(mdata, k, mutual)
  list(graph = igraph::graph_from_edgelist(t(res$edges), directed = !mutual),
       distances = sqrt(res$euc_dists))
}
