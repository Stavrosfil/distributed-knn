#include <iostream>
#include <stdio.h>

class Point {

  public:
    int index;
    double* coords;
    int d;

    Point() : index(-1), coords(nullptr), d(0){};
    Point(int index, double* coords, int d) : index(index), coords(coords), d(d) {}
};

class heapItem {
  public:
    double dist;
    const Point* p;
    heapItem(double dist, const Point* p) : p(p), dist(dist) {}
};

class Node {

  public:
    int vpIndex;
    int parentIndex;
    int leftIndex;
    int rightIndex;
    Point* leafPoints;
    int leafPointsLen;

    double mu;

    Node() : vpIndex(-1), parentIndex(-1), leftIndex(-1), rightIndex(-1), mu(-1.), leafPoints(nullptr), leafPointsLen(0)
    {
    }
};