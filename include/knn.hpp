#ifndef __KNN_H__
#define __KNN_H__

#include <iostream>
#include <stdio.h>
#include <cmath>
#include <limits>
#include <cblas.h>

#define D_MAX std::numeric_limits<double>::max()

typedef struct knnresult {
    int** nidx;     //!< Indices (0-based) of nearest neighbors [m-by-k]
    double** ndist; //!< Distance of nearest neighbors          [m-by-k]
    int m;          //!< Number of query points                 [scalar]
    int k;          //!< Number of nearest neighbors            [scalar]
} knnresult;

void swap(double* a, double* b) {
    double temp = *a;
    *a          = *b;
    *b          = temp;
}

/**
 * @brief Add an element to knn arrays and sort ndist[y][:].
 *
 * @param dist  Distance between point x (X[j]) and y (Y[i])
 * @param x     Index j of X[j]
 * @param ndist kNN distances of Y
 * @param nidx  kNN indices of Y
 * @param k     Number of nearest neighbors
 */
void addTokNN(double dist, int xi, double* ndist, int* nidx, int k) {

    int i    = k - 1;
    ndist[i] = dist;

    while (i > 0 && ndist[i - 1] > ndist[i]) {
        swap(&ndist[i - 1], &ndist[i]);
        i--;
    }
    nidx[i] = xi;
}

void euclideanDistance(double* X, double* Y, double* D, int n, int m, int d) {

    double* XH = new double[n];
    double* YH = new double[m];

    // XH = sum(X.^2,2)
    for (int i = 0; i < n; i++) {
        XH[i] = 0;
        for (int j = 0; j < d; j++) {
            XH[i] += X[i * d + j] * X[i * d + j];
        }
    }

    // YH = sum(Y.^2,2)
    for (int i = 0; i < m; i++) {
        YH[i] = 0;
        for (int j = 0; j < d; j++) {
            YH[i] += Y[i * d + j] * Y[i * d + j];
        }
    }

    // D = 2*X*Y.'
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, m, d, -2, X, d, Y, d, 0, D, m);

    // D = sqrt(sum(X.^2,2) - 2*X*Y.' + sum(Y.^2,2).')
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            D[i * m + j] = sqrt(D[i * m + j] + XH[i] + YH[j]);
        }
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
knnresult kNN(double* X, double* Y, int offset, int n, int m, int d, int k) {

    int** nidx     = new int*[m];
    double** ndist = new double*[m];
    double* D      = new double[n * m];

    for (int i = 0; i < m; i++) {
        nidx[i]  = new int[k];
        ndist[i] = new double[k];
    }

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < k; j++) {
            ndist[i][j] = D_MAX;
        }
    }

    euclideanDistance(X, Y, D, n, m, d);

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double dist = D[j * m + i];
            if (dist == 0)
                continue;
            if (dist < ndist[i][k - 1])
                addTokNN(dist, j + offset, ndist[i], nidx[i], k);
        }
    }

    struct knnresult res = {nidx, ndist, m, n};

    return res;
}

#endif // __KNN_H__