#include <iostream>
#include <stdlib.h>

#include "knn.hpp"
#include "utils.hpp"
#include "distributed.hpp"
#include "vpt_v3.hpp"


int main() {

    // const int n = 5;
    // const int d = 1;
    // const int k = 5;

    // // double X[n * d] = {0, 1, 2, 3, 15, 12, 15, 11, 30, 30};
    // // double X[n * d] = {0, 1, 1, 0, 0, 0, 1, 1, 2, 0};

    // double X[n * d] = {0, 10, -10, 20, 3};

    // struct knnresult ans = mpi::distrAllkNN(X, n, d, k);

  double data[] = {14, 2, 50, 11, 8, 7};

  Node root = Node(0, 6, data, -1);

  // std::cout << "mu = " << root.mu << std::endl;

  VPT t = VPT(root);
  t.createVPT();

  // for(Node n : t.tree) {
  //       std::cout << n.index << '\n';
  // }

  std::cout << "Root:\t\t";
  prt::node(t.tree[0].data, t.tree[0].len);
  std::cout << "Left child:\t";
  prt::node(t.tree[1].data, t.tree[1].len);
  std::cout << "Right child:\t";
  prt::node(t.tree[2].data, t.tree[2].len);

  return 0;
}