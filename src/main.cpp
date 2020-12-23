#include <iostream>
#include <stdlib.h>

#include "knn.hpp"
#include "utils.hpp"
#include "distributed.hpp"
#include "vpt_v2.hpp"


int main() {

    // const int n = 5;
    // const int d = 1;
    // const int k = 5;

    // // double X[n * d] = {0, 1, 2, 3, 15, 12, 15, 11, 30, 30};
    // // double X[n * d] = {0, 1, 1, 0, 0, 0, 1, 1, 2, 0};

    // double X[n * d] = {0, 10, -10, 20, 3};

    // struct knnresult ans = mpi::distrAllkNN(X, n, d, k);

  auto points = std::vector<std::vector<double>> {
    {0, 0, 1},
    {1, 1, 1},
    {2, 0, 0},
    {-1, -1, 0},
    {10, 0, 5}
  };

  vpt::VpTree t1(points); // create a tree

  std::vector<double> distances;
  std::vector<int> indices;
  std::tie(distances, indices) = t1.getNearestNeighbors({ 10, 0, 5 }, 3); // find 3 neighbors closest to the given point

  std::cout << "\nindices" << "\t"; 
  std::cout << "distances" << "\n";
  for (int i = 0; i < 3; i++) {
    std::cout << indices[i] << "\t"; 
    std::cout << distances[i] << "\n"; 
  }
  std::cout << "\n";

  return 0;
}