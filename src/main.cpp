#include <iostream>
#include <stdlib.h>
#include <iterator>
#include <algorithm>

#include "knn.hpp"
#include "utils.hpp"
#include "distributed.hpp"
#include "vpt.hpp"

int main()
{

    std::cout << std::endl;
    // const int n = 5;
    // const int d = 1;
    // const int k = 5;

    // // double X[n * d] = {0, 1, 2, 3, 15, 12, 15, 11, 30, 30};
    // // double X[n * d] = {0, 1, 1, 0, 0, 0, 1, 1, 2, 0};

    // double X[n * d] = {0, 10, -10, 20, 3};

    // struct knnresult ans = mpi::distrAllkNN(X, n, d, k);

    // double data[]      = {100, 80, 70, 60, 40, 35, 200, 500};
    double data[] = {14, 2, 50, 8, 11, 7, 19, 40};
    int k         = 3;
    int d         = 1;
    int len       = 8;
    int m         = 1;

    knnresult ans = knnresult();
    ans.m         = m;
    ans.k         = k;
    ans.nidx      = new int[ans.m * ans.k];
    ans.ndist     = new double[ans.m * ans.k];

    std::fill_n(ans.nidx, ans.m * ans.k, D_MAX);
    std::fill_n(ans.ndist, ans.m * ans.k, -1);

    std::vector<Point> X;

    for (int i = 0; i < len; i += d) {

        double* coords = new double[d];

        for (int j = 0; j < d; j++)
            coords[j] = data[i + j];

        X.push_back(Point(i, coords, d));
    }

    /* ------------------------------ Construct VPT ----------------------------- */

    VPT vpt(X);
    vpt.buildTree(0, X.size());

    int cnt = 0;
    for (auto p : vpt._nodes) {
        std::cout << cnt++ << ": " << X[p.vpIndex].coords[0] << "\t"
                  << "mu: " << p.mu << "\t(" << p.leftIndex << ", " << p.rightIndex << ", " << p.parentIndex << ")"
                  << std::endl;
        if (p.leafPointsLen) {
            for (int j = 0; j < p.leafPointsLen; j++) {
                std::cout << p.leafPoints[j].coords[0] << std::endl;
            }
        }
    }
    std::cout << std::endl;

    Point p = Point(-1, new double[1]{8}, 1);

    /* --------------------------------- VPT kNN -------------------------------- */

    vpt.kNN(p, ans);

    std::cout << "\nkNN ndist:\t";
    prt::rowMajor(ans.ndist, 1, ans.k);
    std::cout << "kNN nidx:\t";
    prt::rowMajor(ans.nidx, 1, ans.k);

    return 0;
}
