#ifdef ENABLE_DOCTEST_IN_LIBRARY
#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"
#endif

#include <iostream>
#include <stdlib.h>

#include "exampleConfig.h"
#include "example.h"

#include "knn.h"

int main() {

    const int n = 2;
    const int m = 3;
    const int d = 2;
    const int k = 2;

    double X[n * d] = {3, 3, 5, 1};
    double Y[m * d] = {0, 1, 2, 5, 4, 1};
    double D[n * m];

    struct knnresult res = kNN(X, Y, D, n, m, d, k);

    std::cout << "kNN distances: " << std::endl;

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < k; j++)
            std::cout << res.ndist[i][j] << " ";
        std::cout << std::endl;
    }

    std::cout << std::endl << "kNN indices: " << std::endl;

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < k; j++)
            std::cout << res.nidx[i][j] << " ";
        std::cout << std::endl;
    }

    return 0;
}