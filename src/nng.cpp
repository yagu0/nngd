#include <iostream>
#include <Rcpp.h>
#include <RcppEigen.h>
#include "knncpp.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
List findNeighbors(NumericMatrix data, int k, bool mutual) {
  int n = data.nrow(),
      d = data.ncol();

  Eigen::MatrixXd dataPoints(d, n);
  for (int row = 0; row < n; ++row) {
    for (int col = 0; col < d; ++col) {
      dataPoints(col,row) = data(row,col); //dataPoints: by columns
    }
  }

  knncpp::KDTreeMinkowskiX<double, knncpp::EuclideanDistance<double>> kdtree(dataPoints);
  kdtree.setBucketSize(16);
  kdtree.setSorted(false);
  kdtree.setTakeRoot(false);
  kdtree.setMaxDistance(0);
  kdtree.setThreads(0);
  kdtree.build();

  knncpp::Matrixi indices;
  Eigen::MatrixXd distances;
  // k+1 because i is always a neighbor of i (to discard)
  kdtree.query(dataPoints, k+1, indices, distances);

  NumericVector res_edges(0);
  NumericVector res_dists(0);
  for (int i = 0; i <= k; ++i) {
    for (int j = 0; j < n; ++j) {
      if (indices(i,j) == j)
        continue;
      bool addRow = false;
      if (!mutual)
        addRow = true;
      else if (mutual && j < indices(i,j)) {
        int l = 0;
        for (; l <= k; ++l) {
          if (indices(l,indices(i,j)) == j)
            break;
        }
        if (l <= k)
          addRow = true;
      }
      if (addRow) {
        // R indices from 1 to n:
        res_edges.push_back(j+1);
        res_edges.push_back(indices(i,j)+1);
        res_dists.push_back(distances(i,j));
      }
    }
  }

  res_edges.attr("dim") = Dimension(2, res_edges.length() / 2);
  List L = List::create(Named("edges") = res_edges,
                        Named("euc_dists") = res_dists);
  return L;
}
