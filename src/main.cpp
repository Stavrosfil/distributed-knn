#include <iostream>
#include <stdlib.h>

#include "knn.hpp"
#include "utils.hpp"
#include "distributed.hpp"

int main() {

    const int n = 5;
    const int d = 1;
    const int k = 5;

    // double X[n * d] = {0, 1, 2, 3, 15, 12, 15, 11, 30, 30};
    // double X[n * d] = {0, 1, 1, 0, 0, 0, 1, 1, 2, 0};

    double X[n * d] = {0, 10, -10, 20, 3};

    struct knnresult ans = mpi::distrAllkNN(X, n, d, k);

    return 0;
}
