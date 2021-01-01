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

struct less_than_key {
    inline bool operator()(const std::pair<double, int> a, const std::pair<double, int> b) {
        return (a.first < b.first);
    }
};

void kNN(knnresult& res, double* X, double* Y, int displacement, int n, int m, int d, int k) {

    double* D = new double[n * m];

    util::computeEuclideanDistance(X, Y, D, n, m, d);

    for (int i = 0; i < m; i++) {

        std::vector<std::pair<double, int>> _D(n + k);

        for (int j = 0; j < n; j++) {
            _D[j].first = D[j * m + i];
            _D[j].second = j + displacement;
        }

        for (int j = 0; j < k; j++) {
            _D[n + j].first = res.ndist[i * k + j];
            _D[n + j].second = res.nidx[i * k + j];
        }

        std::partial_sort(_D.begin(), _D.begin() + k, _D.end(), less_than_key());

        for (int j = 0; j < k; j++) {
            res.ndist[i * k + j] = _D[j].first;
            res.nidx[i * k + j] = _D[j].second;
        }
    }
}

#endif // __KNN_H__