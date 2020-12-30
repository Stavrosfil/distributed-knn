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

    Point(){};
    Point(int index, double* coords, int d) : index(index), coords(coords), d(d) {}
};

class Node {

  public:
    int vpIndex;
    int leftIndex;
    int rightIndex;
    Point* points;

    double mu;

    Node() : vpIndex(-1), leftIndex(-1), rightIndex(-1), mu(0.), points(0) {}
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

    VPT(std::vector<Point> points) : _points(points) {}

    // Partition and construct tree from _points[lo:hi]
    // lo and hi are indices of corpus array
    int buildTree(int lo, int hi) {

        if (lo == hi)
            return -1;

        Node node;
        node.vpIndex = lo;

        if (hi - lo > 1) {
            int median = (hi + lo) % 2 == 0 ? (hi + lo) / 2 : (hi + lo + 1) / 2;

            node.mu         = computeMu(lo, median, hi);
            node.leftIndex  = buildTree(lo + 1, median);
            node.rightIndex = buildTree(median, hi);
        }

        _nodes.push_back(node);
        return _nodes.size() - 1;
    }

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