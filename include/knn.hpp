#ifndef __KNN_H__
#define __KNN_H__

#include <iostream>
#include <stdio.h>
#include <cmath>
#include <cblas.h>
#include <algorithm>
#include <vector>

#define D_MAX std::numeric_limits<double>::max()

#include "utils.hpp"

typedef struct knnresult {
    int* nidx;     //!< Indices (0-based) of nearest neighbors [m-by-k]
    double* ndist; //!< Distance of nearest neighbors          [m-by-k]
    int m;         //!< Number of query points                 [scalar]
    int k;         //!< Number of nearest neighbors            [scalar]
} knnresult;

void swap(double* a, double* b) {
    double temp = *a;
    *a          = *b;
    *b          = temp;
}

void swap(int* a, int* b) {
    int temp = *a;
    *a       = *b;
    *b       = temp;
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

    // D = -2*X*Y.'
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, m, d, -2, X, d, Y, d, 0, D, m);

    // D = sqrt(sum(X.^2,2) - 2*X*Y.' + sum(Y.^2,2).')
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            D[i * m + j] = sqrt(D[i * m + j] + XH[i] + YH[j]);
        }
    }
}

struct less_than_key {
    inline bool operator()(const std::pair<double, int> a, const std::pair<double, int> b) {
        return (a.first < b.first);
    }
};

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
void kNN(knnresult res, double* X, double* Y, int displacement, int n, int m, int d, int k) {

    // std::cout << "X sub-matrix in kNN:\n";
    // prt::rowMajor(Y, n, d);
    // std::cout << "\n";

    double* D = new double[n * m];

    euclideanDistance(X, Y, D, n, m, d);

    for (int i = 0; i < m; i++) {

        std::vector<std::pair<double, int>> _D;
        _D.reserve(n + k);

        for (int j = 0; j < n; j++) {
            _D[j].first  = D[j * m + i];
            _D[j].second = j + displacement;
        }

        for (int j = 0; j < k; j++) {
            _D[n + j].first  = res.ndist[i * k + j];
            _D[n + j].second = res.nidx[i * k + j];
        }

        std::partial_sort(_D.begin(), _D.begin() + k, _D.begin() + n + k, less_than_key());
        // std::nth_element(_D.begin() + (col), _D.begin() + (col + k - 1), _D.begin() + (col + m), less_than_key());

        for (int j = 0; j < k; j++) {
            res.ndist[i * k + j] = _D[j].first;
            res.nidx[i * k + j]  = _D[j].second;
        }
    }
}

#endif // __KNN_H__