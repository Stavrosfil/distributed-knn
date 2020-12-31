#include <iostream>
#include <stdio.h>


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
    int parentIndex;
    int leftIndex;
    int rightIndex;
    Point* points;
    int points_len;

    double mu;

    Node() : vpIndex(-1), parentIndex(-1), leftIndex(-1), rightIndex(-1), mu(-1.), points(nullptr), points_len(0) {}
};