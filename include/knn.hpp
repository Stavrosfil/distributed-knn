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
typedef std::priority_queue<heapItem, std::vector<heapItem>, comp::heapDist> pointHeap;

#define Y_BLOCKS 1000

/* -------------------------------- kNN function -------------------------------- */

// Finds for each point in a query set Y the k nearest neighbors in the corpus set X
// X[n * d]
// Y[m * d]
// res.nidx[m * k]
// res.dist[m * k]

void kNN(knnresult& res, double* X, double* Y, int displacement, int n, int m, int d, int k)
{
    std::vector<int> chunk_size(Y_BLOCKS);
    std::vector<int> displs(Y_BLOCKS);

    util::computeChunksDisplacements(chunk_size.data(), displs.data(), Y_BLOCKS, m, d);

    int cnt = 0;
    for (auto displ : displs) {

        displ /= d;
        int y_len = chunk_size[cnt++] / d;

        std::vector<double> D(n * y_len, 0);
        util::computeEuclideanDistance(X, Y + displ * d, D.data(), n, y_len, d);

        for (int i = displ; i < displ + y_len; i++) {

            std::vector<std::pair<double, int>> _D(n + k);

            for (int j = 0; j < n; j++) {
                _D[j].first  = D[j * y_len + i - displ];
                _D[j].second = j + displacement;
            }

            for (int j = 0; j < k; j++) {
                _D[n + j].first  = res.ndist[i * k + j];
                _D[n + j].second = res.nidx[i * k + j];
            }

            std::partial_sort(_D.begin(), _D.begin() + k, _D.end(), comp::lessThanKey());

            for (int j = 0; j < k; j++) {
                res.ndist[i * k + j] = _D[j].first;
                res.nidx[i * k + j]  = _D[j].second;
            }
        }
    }
}

/* ---------------------------- update kNN function ----------------------------- */

// Used to update kNN in leafs of VPT
// TODO optimize (use euclideanDistance)

void updateKNN(pointHeap& heap, Point& queryPoint, Point& corpusPoint, int k)
{
    double dist = util::distance(queryPoint, corpusPoint);
    if (heap.size() == k) {
        if (dist < heap.top().dist) {
            heap.pop();
            heap.push(heapItem(dist, corpusPoint.index));
        }
    }
    else {
        heap.push(heapItem(dist, corpusPoint.index));
    }
}

#endif // __KNN_H__