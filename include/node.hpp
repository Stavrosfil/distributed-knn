#include <iostream>
#include <stdio.h>


class Point {

  public:
    uint32_t index;
    double* coords;
    int d;

    Point() : index(-1), coords(nullptr), d(0){};
    Point(int index, double* coords, int d) : index(index), coords(coords), d(d) {}
};

class Node {

  public:
    int vpIndex;
    int parentIndex;
    int leftIndex;
    int rightIndex;
    Point* points;
    int points_len;

    double mu;

    Node() : vpIndex(-1), parentIndex(-1), leftIndex(-1), rightIndex(-1), mu(-1.), points(nullptr), points_len(0) {}
};