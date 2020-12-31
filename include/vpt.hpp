#include <iostream>
#include <stdio.h>
#include <cmath>
#include <array>
#include <vector>

#include "knn.hpp"

class Point {

  public:
    uint32_t index;
    uint8_t d;
    double* coords;

    Point() : index(-1), d(0), coords(nullptr){};
    Point(int index, double* coords, int d) : index(index), coords(coords), d(d) {}
};

class Node {

  public:
    int vpIndex;
    int leftIndex;
    int rightIndex;
    Point* points;
    int points_len;

    double mu;

    Node() : vpIndex(-1), leftIndex(-1), rightIndex(-1), mu(-1.), points(nullptr), points_len(0) {}
};

double distance(Point& p1, Point& p2) {
    double* D = new double[1];
    util::computeEuclideanDistance(p1.coords, p2.coords, D, 1, 1, p1.d);
    return D[0];
}

class VPT {

  public:
    std::vector<Point> _points = {};
    std::vector<Node> _nodes   = {};
    int b                      = 0;

    VPT(std::vector<Point> points) : _points(points) {}

    // Partition and construct tree from _points[lo:hi]
    // lo and hi are indices of corpus array
    int buildTree(int lo, int hi) {

        if (lo == hi)
            return -1;

        Node node;
        node.vpIndex = lo;

        int n = hi - (lo + 1);
        if (n > 2 * b + 1 || b == 0) {
            int median      = (hi + lo) % 2 == 0 ? (hi + lo) / 2 : (hi + lo + 1) / 2;
            node.mu         = computeMu(lo + 1, median, hi);
            node.leftIndex  = buildTree(lo + 1, median);
            node.rightIndex = buildTree(median, hi);
        } else {
            node.points     = new Point[n];
            node.points_len = n;
            for (int i = 0; i < n; i++) {
                node.points[i] = _points[lo + 1 + i];
            }
        }

        _nodes.push_back(node);
        return _nodes.size() - 1;
    }

  private:
    double computeMu(int lo, int median, int hi) {
        std::nth_element(_points.begin() + lo + 1, _points.begin() + median, _points.begin() + hi,
                         distanceFromVP(_points[lo]));
        return distance(_points[lo], _points[median]);
    }

    struct distanceFromVP {
        Point& p;
        distanceFromVP(Point& p) : p(p) {}
        bool operator()(Point& p1, Point& p2) {
            return distance(p, p1) < distance(p, p2);
        }
    };
};