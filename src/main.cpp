#include <iostream>
#include <stdlib.h>

#include "knn.hpp"
#include "utils.hpp"
#include "distributed.hpp"
#include "vpt.hpp"


int main() {

    std::cout << std::endl;
    // const int n = 5;
    // const int d = 1;
    // const int k = 5;

    // // double X[n * d] = {0, 1, 2, 3, 15, 12, 15, 11, 30, 30};
    // // double X[n * d] = {0, 1, 1, 0, 0, 0, 1, 1, 2, 0};

    // double X[n * d] = {0, 10, -10, 20, 3};

    // struct knnresult ans = mpi::distrAllkNN(X, n, d, k);

  double data[] = {14, 2, 50, 8, 11, 7, 19, 40};
  int len = 8;
  int k = 2;
  int d = 1;
  double queryPoints = 19;
  int m = 1;
  Node root = Node(0, len, data, -1);

  // std::cout << "mu = " << root.mu << std::endl;

  VPT t = VPT(root);
  t.createVPT();

  knnresult _ans = knnresult();
  _ans.m         = m;
  _ans.k         = k;
  _ans.nidx      = new int[_ans.m * _ans.k]();
  _ans.ndist     = new double[_ans.m * _ans.k]();

  for (int i = 0; i < _ans.m; i++) {
      for (int j = 0; j < k; j++) {
          _ans.ndist[i * k + j] = D_MAX;
          _ans.nidx[i * k + j]  = -1;
      }
  }
  
  t.vptKnn(_ans, &queryPoints, 0, _ans.m, d, k);

  std::cout << "\nkNN ndist: ";
  prt::rowMajor(_ans.ndist, 1, k);
  std::cout << "\n";

  for(Node n : t.tree) {
    std::cout << "Node\t" << n.index << ":\t\t";
    prt::node(n.data, n.len);
    // std::cout << "parent = " << n.parentIndex << "\tleft = " << n.leftIndex << "\tright = " << n.rightIndex << "\tmu = " << n.mu << std::endl;
  }
  std::cout << std::endl;
  
  return 0;
}