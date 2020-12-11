#ifdef ENABLE_DOCTEST_IN_LIBRARY
#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"
#endif

#include <iostream>
#include <stdlib.h>

#include "exampleConfig.h"
#include "example.h"

#include "knn.h"

/*
 * Simple main program that demontrates how access
 * CMake definitions (here the version number) from source code.
 */
int main() {

    int n = 7;
    int m = 3;
    int d = 1;
    int k = 3;

    // double **X = new double *[n];
    // double **Y = new double *[m];

    // for (int i = 0; i < n; i++) {
    //     X[i] = new double[d];
    // }

    // for (int i = 0; i < m; i++) {
    //     Y[i] = new double[d];
    // }

    double X[][1] = {{1.0}, {2.0}, {4.0}, {15.0}, {16.0}, {18.0}, {40.0}};
    double Y[][1] = {{0.0}, {17.0}, {30.0}};

    struct knnresult res = kNN(*X, *Y, n, m, d, k);

    std::cout << "kNN distances: " << std::endl;

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < k; j++)
            std::cout << res.ndist[i][j] << " ";
        std::cout << std::endl;
    }

    std::cout << std::endl << "kNN indeces: " << std::endl;

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < k; j++)
            std::cout << res.nidx[i][j] << " ";
        std::cout << std::endl;
    }

    return 0;
}