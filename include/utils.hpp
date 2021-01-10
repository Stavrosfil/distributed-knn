#ifndef __UTILS_H__
#define __UTILS_H__

#include <iostream>
#include <stdlib.h>
#include "node.hpp"

struct knnresult {
    int* nidx;     //!< Indices (0-based) of nearest neighbors [m-by-k]
    double* ndist; //!< Distance of nearest neighbors          [m-by-k]
    int m;         //!< Number of query points                 [scalar]
    int k;         //!< Number of nearest neighbors            [scalar]
};

namespace util {

class Timer {
  public:
    Timer(bool print) : print(print) {}

    void start(std::string operation_desc)
    {
        _operation_desc = operation_desc;
        t1              = std::chrono::high_resolution_clock::now();
    }

    void stop()
    {
        t2       = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        if (print)
            std::cout << "\n" << _operation_desc << " time: " << duration / 1e3 << "ms\n" << std::endl;
    }

  private:
    double duration;
    bool print;
    std::string _operation_desc;
    std::chrono::high_resolution_clock::time_point t1, t2;
};

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

    std::vector<double> XH(n, 0);
    std::vector<double> YH(m, 0);

    // XH = sum(X.^2,2)
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < d; j++) {
            XH[i] += X[i * d + j] * X[i * d + j];
        }
    }

    // YH = sum(Y.^2,2)
    for (int i = 0; i < m; i++) {
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

void readToRowMajorVector(std::vector<double>& result, int n, int d, std::string fileName)
{
    std::cout << fileName << std::endl;
    std::string line;
    std::string filePath = "./dataset/" + fileName;
    std::ifstream myfile(filePath);

    std::ifstream input(filePath);

    std::string s;

    std::getline(input, s);
    std::istringstream iss(s);
    std::cout << s << std::endl;

    for (int i = 0; i < n; i++) {
        std::getline(input, s);
        std::istringstream iss(s);

        std::string num;
        int j = 0;
        while (std::getline(iss, num, ',')) {
            result[i * d + j++] = std::stof(num);
        }
    }
}

void read(std::vector<double>& result, int n, int d, std::string fileName)
{
    std::cout << fileName << std::endl;
    std::string line;
    std::string filePath = "./dataset/" + fileName;
    std::ifstream myfile(filePath);

    std::ifstream input(filePath);

    std::string s;

    std::getline(input, s);
    std::istringstream iss(s);
    std::cout << s << std::endl;

    for (int i = 0; i < n; i++) {
        std::getline(input, s);
        std::istringstream iss(s);

        std::string num;
        int j     = 0;
        int count = 0;
        while (std::getline(iss, num, ' ')) {
            if (count++)
                result[i * d + j++] = std::stof(num);
        }
    }
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

void vector(std::vector<int> v)
{
    for (int i = 0; i < v.size(); i++) {
        std::cout << v.data()[i] << "\t";
    }
    std::cout << std::endl;
}

void vector(std::vector<double> v)
{
    for (int i = 0; i < v.size(); i++) {
        std::cout << v.data()[i] << "\t";
    }
    std::cout << std::endl;
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
    std::cout << p.index << "( ";
    for (int i = 0; i < p.d; i++)
        std::cout << p.coords[i] << " ";
    std::cout << ")\t";
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

void kNN(knnresult res)
{
    std::cout << "\nkNN results:\n\n";
    for (int i = 0; i < res.m; i++) {
        std::cout << "nidx[" << i << "]"
                  << ":\t";
        for (int j = 0; j < res.k; j++) {
            std::cout << res.nidx[i * res.k + j] << "\t";
        }
        std::cout << "\nndist[" << i << "]"
                  << ":\t";
        for (int j = 0; j < res.k; j++) {
            std::cout << res.ndist[i * res.k + j] << "\t";
        }
        std::cout << std::endl << std::endl;
    }
}

} // namespace prt

namespace conv {

void serVector(std::vector<Point>& points, std::vector<int>& _indices, std::vector<double>& _coords, int len)
{
    int _d = points[0].d;

    for (int i = 0; i < len; i++) {
        for (int j = 0; j < _d; j++) {
            _coords[i * _d + j] = points[i].coords[j];
        }
        _indices[i] = points[i].index;
    }
}

void recVector(std::vector<Point>& points, std::vector<int>& _indices, std::vector<double>& _coords, int d, int len)
{
    for (int i = 0; i < len; i++) {
        points[i].d      = d;
        points[i].index  = _indices[i];
        points[i].coords = new double[d];
        for (int j = 0; j < d; j++) {
            points[i].coords[j] = _coords[i * d + j];
        }
    }
}

} // namespace conv

#endif // __UTILS_H__