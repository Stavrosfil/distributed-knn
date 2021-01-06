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
    int index;

    heapItem(double dist, int index) : dist(dist), index(index) {}
};

class Node {

  public:
    int vpIndex;
    double mu;

    Node* left;
    Node* right;

    int leafPointsIndex;
    int leafPointsLen;

    Node() : vpIndex(-1), left(nullptr), right(nullptr), mu(-1.), leafPointsIndex(-1), leafPointsLen(0) {}
};