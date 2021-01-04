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

double computeEucledianNorm(const double* d1, const double* d2, int d)
{
    // double* dist = new double[d];
    double dist = 0;
    for (int i = 0; i < d; i++) {
        dist += (d1[i] - d2[i]) * (d1[i] - d2[i]);
    }
    // double ans = cblas_dnrm2(d, dist, 1);
    // delete dist;
    return sqrt(dist);
}

double distance(const Point& p1, const Point& p2)
{
    return util::computeEucledianNorm(p1.coords, p2.coords, p1.d);
}

} // namespace util

namespace comp {

struct heapDist {
    bool operator()(const heapItem& lhs, const heapItem& rhs)
    {
        return lhs.dist < rhs.dist;
    }
};

struct lessThanKey {
    inline bool operator()(const std::pair<double, int>& a, const std::pair<double, int>& b)
    {
        return (a.first < b.first);
    }
};

struct distanceFromVP {
    const Point& p;
    distanceFromVP(const Point& p) : p(p) {}
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

void point(Point& p)
{
    std::cout << "( ";
    for (int i = 0; i < p.d; i++)
        std::cout << p.coords[i] << " ";
    std::cout << ")\t";
    // std::cout << std::endl << std::endl;
}

void tree(Node* root, std::vector<Point>& points)
{
    std::vector<std::vector<int>> levels;
    if (root == NULL) {
        return;
    }
    std::vector<int> level;
    std::queue<Node*> q;
    q.push(root);

    while (!q.empty()) {
        int curr_size = q.size();
        for (int i = 0; i < curr_size; ++i) {
            Node* front = q.front();
            q.pop();
            prt::point(points[front->vpIndex]);
            std::cout << " (mu = " << front->mu << ") ->\t";
            level.push_back(front->vpIndex);
            if (front->left) {
                q.push(front->left);
                prt::point(points[front->left->vpIndex]);
            }
            if (front->right) {
                q.push(front->right);
                prt::point(points[front->right->vpIndex]);
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
        levels.push_back(level);
        level.clear();
    }
}

} // namespace prt

namespace conv {

void serVector(std::vector<Point>& points, int* _indices, double* _coords)
{
    int _len = points.size();
    int _d = points[0].d;

    for(int i = 0; i < _len; i++) {
        for(int j = 0; j < _d; j++) {
            _coords[i * _d + j] = points[i].coords[j];
        }
        _indices[i] = points[i].index;
    }
}

void recVector(std::vector<Point>& points, int* _indices, double* _coords, int len, int d)
{
    for (int i = 0; i < len; i++) {
        points[i].coords = new double[d];
        points[i].d = d;
        for (int j = 0; j < d; j++) {
            points[i].coords[j] = _coords[i * d + j];
        }
        points[i].index = _indices[i];
    }
}

}

#endif // __UTILS_H__