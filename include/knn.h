#ifndef KNN
#define KNN
#include <iostream>
#include <stdio.h>
#include <cmath>
#include <limits>

#define D_MAX std::numeric_limits<double>::max()

/**
 * @brief kNN result struct
 *
 */
typedef struct knnresult {
    int **nidx;     //!< Indices (0-based) of nearest neighbors [m-by-k]
    double **ndist; //!< Distance of nearest neighbors          [m-by-k]
    int m;          //!< Number of query points                 [scalar]
    int k;          //!< Number of nearest neighbors            [scalar]
} knnresult;

/**
 * @brief Add an element to knn arrays and sort ndist[y][:].
 *
 * @param dist  Distance between point x (X[j]) and y (Y[i])
 * @param x     X[j]
 * @param ndist kNN distances of Y
 * @param nidx  kNN indeces of Y
 * @param k     Number of nearest neighbors to search for
 */
void addTokNN(double dist, int x, double *ndist, int *nidx, int k) {

    // TODO: Binary search for insertion

    for (int i = k - 1; i > 0; i--) {
        if (ndist[i - 1] < dist) {
            for (int j = k - 1; j > i; j--) {
                ndist[j] = ndist[j - 1];
                nidx[j]  = nidx[j - 1];
            }
            ndist[i] = dist;
            nidx[i]  = x;
            break;
        }
    }
    if (ndist[0] >= dist) {
        for (int j = k - 1; j > 0; j--) {
            ndist[j] = ndist[j - 1];
            nidx[j]  = nidx[j - 1];
        }
        ndist[0] = dist;
        nidx[0]  = x;
    }
}

/**
 * @brief Compute k nearest neighbors of each point in X [n-by-d]
 *
 * \param  X      Corpus data points              [n-by-d]
 * \param  Y      Query data points               [m-by-d]
 * \param  n      Number of corpus points         [scalar]
 * \param  m      Number of query points          [scalar]
 * \param  d      Number of dimensions            [scalar]
 * \param  k      Number of neighbors             [scalar]
 *
 * \return  The kNN result
 */
knnresult kNN(double *X, double *Y, int n, int m, int d, int k) {

    int **nidx     = new int *[m];
    double **ndist = new double *[m];

    for (int i = 0; i < m; i++) {
        nidx[i]  = new int[k];
        ndist[i] = new double[k];
    }

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < k; j++) {
            ndist[i][j] = D_MAX;
        }
    }

    // TODO: Use BLAS for distance function for n-dimensions.
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double dist = abs(Y[i * d + 0] - X[j * d + 0]);
            if (dist == 0)
                continue;
            if (dist < ndist[i][k - 1]) 
                addTokNN(dist, j, ndist[i], nidx[i], k);
            
        }
    }

    struct knnresult res = {nidx, ndist, m, n};

    return res;
}

#endif