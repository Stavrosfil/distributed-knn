#ifndef __UTILS_H__
#define __UTILS_H__

#include <iostream>
#include <stdlib.h>
#include "node.hpp"

namespace util {

    void computeChunksDisplacements(int* cnt, int* displs, int processes, int row_s, int col_s)
    {
        int rem = row_s % processes; // elements remaining after division among processes
        int sum = 0;                 // Sum of counts. Used to calculate displacements

        for (int i = 0; i < processes; i++) {
            cnt[i] = (row_s / processes) * col_s;
            if (rem > 0) {
                cnt[i] += col_s;
                rem--;
            }
            displs[i] = sum;
            sum += cnt[i];
        }
    }

    void computeEuclideanDistance(double* X, double* Y, double* D, int n, int m, int d)
    {

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

    double distance(const Point& p1, const Point& p2)
    {
        double* D = new double[1];
        util::computeEuclideanDistance(p1.coords, p2.coords, D, 1, 1, p1.d);
        return D[0];
    }

} // namespace util

namespace comp {

    struct heapDist {
        bool operator()(const std::pair<double, Point>& lhs, const std::pair<double, Point>& rhs)
        {
            return lhs.first < rhs.first;
        }
    };

    struct lessThanKey {
        inline bool operator()(const std::pair<double, int> a, const std::pair<double, int> b)
        {
            return (a.first < b.first);
        }
    };

    struct distanceFromVP {
        Point& p;
        distanceFromVP(Point& p) : p(p) {}
        bool operator()(Point& p1, Point& p2)
        {
            return util::distance(p, p1) < util::distance(p, p2);
        }
    };

} // namespace comp

namespace prt {

    void rowMajor(double* a, int n, int m)
    {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++)
                std::cout << a[i * m + j] << "\t";
            std::cout << std::endl;
        }
    }

    void rowMajor(int* a, int n, int m)
    {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++)
                std::cout << a[i * m + j] << "\t";
            std::cout << std::endl;
        }
    }

    void rowMajor(std::pair<double, int>* a, int n, int m)
    {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++)
                std::cout << a[i * m + j].first << "\t";
            std::cout << std::endl;
        }
    }

    void points(std::vector<Point> _points)
    {
        for (int i = 0; i < _points.size(); i++) {
            std::cout << _points[i].index << " -> ( ";
            for (int j = 0; j < _points[0].d; j++)
                std::cout << _points[i].coords[j] << " ";
            std::cout << ")\t";
        }
        std::cout << std::endl << std::endl;
    }

} // namespace prt

#endif // __UTILS_H__