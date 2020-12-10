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

    // int n = 2;
    // int m = 3;
    // int d = 2;
    // int k = 2;
    // int X[n][d] = {{3, 3}, {5, 1}};
    // int Y[m][d] = {{0, 1}, {2, 5}, {4, 1}};

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

    kNN(*X, *Y, n, m, d, k);

    return 0;
}