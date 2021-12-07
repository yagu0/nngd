#' Compute the Euclidian Commute-Time Distances (ECTD) in an undirected graph.
#'
#' Assuming similarity function doesn't depend on x, and undirected graph.
#'
#' @param o Output of \code{nng}.
#' @param similarity function distance --> similarity.
#' @param inf.replace Replace Inf values by large finite number ?
#'
#' @return A distances matrix (n x n)
#'
#' @export
ectd <- function(o, similarity=function(x) 1, inf.replace=TRUE)
{
  if (igraph::is_directed(o$graph))
    stop("Directed graph case unimplemented yet")
  n <- igraph::vcount(o$graph)
  E <- igraph::as_edgelist(o$graph, names=FALSE)
  W <- matrix(0, nrow=n, ncol=n)
  for (i in seq_len(nrow(E))) {
    edge <- c(E[i,1], E[i,2])
    W[edge[1], edge[2]] <- similarity(o$distances[i])
    W[edge[2], edge[1]] <- W[edge[1], edge[2]]
  }
  L <- diag(rowSums(W)) - W
  cc <- igraph::components(o$graph)
  distances <- matrix(Inf, nrow=n, ncol=n)
  for (i in seq_len(cc$no)) {
    indices <- which(cc$membership == i)
    L_loc <- L[indices, indices]
    n_loc <- cc$csize[i]
    if (n_loc >= 2) L_inv = pracma::pinv(L_loc)
    else L_inv = matrix(ifelse(L_loc != 0, 1/L_loc[1], 0))
    Lii = matrix(rep(diag(L_inv), each=n_loc), nrow=n_loc, byrow=TRUE)
    Ljj = matrix(rep(diag(L_inv), each=n_loc), nrow=n_loc)
    # https://math.stackexchange.com/questions/1321305/commute-time-distance-in-a-graph
    distances[indices, indices] <- Lii + Ljj - 2 * L_inv
  }
  if (inf.replace) {
    finiteMax <- max(distances[is.finite(distances)])
    distances[is.infinite(distances)] <- finiteMax + 1
  }
  distances
}
