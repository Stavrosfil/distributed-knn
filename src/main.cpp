#include <iostream>
#include <stdlib.h>

#include "knn.h"

/*
 * Simple main program that demontrates how access
 * CMake definitions (here the version number) from source code.
 */
int main() {

    int n = 2;
    int m = 3;
    int d = 2;
    int k = 2;
    double X[][d] = {{3, 3}, {5, 1}};
    double Y[][d] = {{0, 1}, {2, 5}, {4, 1}};
    double D[n][m];

    struct knnresult res = kNN(*X, *Y, *D, n, m, d, k);

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