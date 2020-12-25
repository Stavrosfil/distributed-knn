#include <iostream>
#include <stdlib.h>

#include "knn.hpp"
#include "utils.hpp"
#include "distributed.hpp"
#include "vpt.hpp"


int main() {

    // const int n = 5;
    // const int d = 1;
    // const int k = 5;

    // // double X[n * d] = {0, 1, 2, 3, 15, 12, 15, 11, 30, 30};
    // // double X[n * d] = {0, 1, 1, 0, 0, 0, 1, 1, 2, 0};

    // double X[n * d] = {0, 10, -10, 20, 3};

    // struct knnresult ans = mpi::distrAllkNN(X, n, d, k);

  double data[] = {14, 2, 50, 8, 11, 7, 19, 40};
  int len = 8;
  Node root = Node(0, len, data, -1);

  // std::cout << "mu = " << root.mu << std::endl;

  VPT t = VPT(root);
  t.createVPT();

  // struct knnresult ans;

  // t.vptKnn(ans, data, data[2], 0, 8, 1, 1, 1);

  for(Node n : t.tree) {
    std::cout << "Node\t" << n.index << ":\t\t";
    // prt::node(n.data, n.len);
    std::cout << "parent = " << n.parentIndex << "\tleft = " << n.leftIndex << "\tright = " << n.rightIndex << "\tmu = " << n.mu << std::endl;
  }

 
  // std::cout <<  t.moveRight(t.tree[0].index) << std::endl;

  // std::cout << t.searchLeaf(14) << std::endl;

  return 0;
}