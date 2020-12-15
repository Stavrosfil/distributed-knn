#include <iostream>
#include <stdlib.h>

#include "knn.hpp"
#include "utils.hpp"
#include "distributed.hpp"

int main() {

    /* ---------------------------------- DATA ---------------------------------- */

    const int n = 5;
    const int m = 3;
    const int d = 1;

    // double X[n * d] = {3, 3, 5, 1};
    // double Y[m * d] = {0, 1, 2, 5, 4, 1};

    double X[n * d] = {0, 10, -10, 20, 30};
    double Y[m * d] = {1, 3, 9};

    double D[m * n];

    const int k = 2;

    /* ----------------------------------- kNN ---------------------------------- */

    struct knnresult res = kNN(X, Y, D, n, m, d, k);

    /* ----------------------------------- MPI ---------------------------------- */

    mpi::test();

    /* --------------------------------- PRINTS --------------------------------- */

    // std::cout << "kNN distances: " << std::endl;
    // prt::twoDim(res.ndist, m, k);

    // std::cout << std::endl << "kNN indices: " << std::endl;
    // prt::twoDim(res.nidx, m, k);

    // std::cout << std::endl << "Distance matrix: " << std::endl;
    // prt::rowMajor(D, n, m);

    return 0;
}
