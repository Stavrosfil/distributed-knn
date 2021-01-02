#ifndef __KNN_H__
#define __KNN_H__

#include <iostream>
#include <stdio.h>
#include <cmath>
#include <cblas.h>
#include <algorithm>
#include <vector>
#include <queue>

#include "utils.hpp"

#define D_MAX std::numeric_limits<double>::max()
typedef std::priority_queue<std::pair<double, Point>, std::vector<std::pair<double, Point>>, comp::heapDist> pointHeap;

struct knnresult {
    int* nidx;     //!< Indices (0-based) of nearest neighbors [m-by-k]
    double* ndist; //!< Distance of nearest neighbors          [m-by-k]
    int m;         //!< Number of query points                 [scalar]
    int k;         //!< Number of nearest neighbors            [scalar]
};

namespace prt {

    void kNN(knnresult res)
    {
        std::cout << "\nkNN results:\n\n";
        for (int i = 0; i < res.m; i++) {
            std::cout << "nidx[" << i << "]" << ":\t";
            for (int j = 0; j < res.k; j++) {
                std::cout << res.nidx[i * res.k + j] << "\t";
            }
            std::cout << "\nndist[" << i << "]" << ":\t";
            for (int j = 0; j < res.k; j++) {
                std::cout << res.ndist[i * res.k + j] << "\t";
            }
            std::cout << std::endl << std::endl;
        }
    }
}

void kNN(knnresult& res, double* X, double* Y, int displacement, int n, int m, int d, int k)
{
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

        std::partial_sort(_D.begin(), _D.begin() + k, _D.end(), comp::lessThanKey());

        for (int j = 0; j < k; j++) {
            res.ndist[i * k + j] = _D[j].first;
            res.nidx[i * k + j] = _D[j].second;
        }
    }
}

void updateKNN(pointHeap& heap, Point& queryPoint, Point& corpusPoint, int k)
{
    double dist = util::distance(queryPoint, corpusPoint);
    if (heap.size() == k) {
        if (dist < heap.top().first) {
            heap.pop();
            heap.push(std::make_pair(dist, corpusPoint));
        }
    }
    else {
        heap.push(std::make_pair(dist, corpusPoint));
    }
}

#endif // __KNN_H__