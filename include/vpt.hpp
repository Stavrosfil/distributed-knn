#include <iostream>
#include <stdio.h>
#include <cmath>
#include <array>
#include <vector>

#include "knn.hpp"

class Point {

  public:
    int index;
    double* coords;

    Point(int d) {
        coords = new double[d];
    }

    Point(double* d) {
        coords = d;
    }

    ~Point() {
        delete coords;
    }
};

class Node {

  public:
    int vpIndex;
    double mu;
    Node* left;
    Node* right;

    Node() : vpIndex(0), mu(0.), left(0), right(0) {}

    ~Node() {
        delete left;
        delete right;
    }
};

class VPT {

  public:
    // double* D                 = {};
    Node* root;
    std::vector<double> D     = {};
    std::vector<Point> points = {};
    int d                     = 1;

    VPT(int d) : d(d) {}

    void computeDistances() {
        util::computeEuclideanDistance(data, &vp, D.data(), len, 1, 1);
    }

    // Partition and construct tree from X[lo:hi]
    Node* buildTree(int lo, int hi) {

        if (lo == hi)
            return NULL;

        Node* node    = new Node();
        node->vpIndex = lo;

        if (hi - lo > 1) {

            int median = (hi + lo) / 2;

            std::nth_element(points.begin() + lo, points.begin() + median, points.begin() + hi + 1,
                             distanceFromVP(node->vpIndex));
        }
    }

    struct distanceFromVP {
        // Point* p;
        // distanceFromVP(int index) {
        //     p = new Point(index);
        // }
        bool operator()(Point* p1, Point* p2) {
            return D[p1->index];
        }
    };

    // double computeMu() {

    //     double* sortedD = new double[len];
    //     memcpy(sortedD, D, len * sizeof(double));
    //     qsort(sortedD, len, sizeof(double), util::compare); // TODO more efficient
    //     mu = sortedD[(len - 1) / 2];
    //     // prt::rowMajor(sortedD, 1, len);
    //     delete[] sortedD;
    //     return mu;
    // }

    void createVPT() {

        int len    = -1;
        int maxIdx = 2 * root.len - 1;
        // std::cout << "Max Idx = " << maxIdx << std::endl;
        for (int i = 0; i < maxIdx; i++) {
            len = tree[i].len;
            if (len > 1) {
                tree[i].leftIndex  = tree.size();
                tree[i].rightIndex = tree.size() + 1;
                createChildren(tree[i], tree[i].len);
                // delete[] tree[i].data;
            }
        }
    }

    void createChildren(Node curNode, int len) {

        int leftLen  = (curNode.len + 1) / 2;
        int rightLen = curNode.len - leftLen;
        Node _left   = Node(tree.size(), leftLen, computeLeftChild(&curNode), curNode.index);
        tree.push_back(_left);
        Node _right = Node(tree.size(), rightLen, computeRightChild(&curNode), curNode.index);
        tree.push_back(_right);
    }

    bool isLeaf(int nodeIdx) {
        if (tree[nodeIdx].leftIndex == -1)
            return true;
        else
            return false;
    }

    int moveUp(int curNodeIndex) {
        return tree[curNodeIndex].parentIndex;
    }

    int moveLeft(int curNodeIndex) {
        return tree[curNodeIndex].leftIndex;
    }

    int moveRight(int curNodeIndex) {
        return tree[curNodeIndex].rightIndex;
    }

    void vptKnn(knnresult res, double* y, int displacement, int m, int d, int k) {

        // for i = 0 -> i < m
        int curNodeIndex = tree[0].index;
        int curNodeLen   = tree[0].len;
        double* tau      = &res.ndist[k - 1];
        // std::cout << *tau << std::endl;

        if (checkInside(tree[curNodeIndex].vp, *y, tree[curNodeIndex].mu, *tau)) {
            searchLeftSubtree(res, curNodeIndex, *y, *tau, displacement, d, k);
        }
        if (checkOutside(tree[curNodeIndex].vp, *y, tree[curNodeIndex].mu, *tau)) {
            searchRightSubtree(res, curNodeIndex, *y, *tau, displacement, d, k);
        }
    }

    void searchLeftSubtree(knnresult res, int curNodeIndex, double y, double tau, int displacement, int d, int k) {

        std::cout << "Node " << curNodeIndex
                  << " inside condition = " << checkInside(tree[curNodeIndex].vp, y, tree[curNodeIndex].mu, tau)
                  << " -> search left subtree\n";
        curNodeIndex = moveLeft(curNodeIndex);
        if (isLeaf(curNodeIndex)) {
            std::cout << "Reached node " << curNodeIndex << std::endl << std::endl;
            kNN(res, tree[curNodeIndex].data, &y, displacement, tree[curNodeIndex].len, 1, d, k);
        } else {
            if (checkInside(tree[curNodeIndex].vp, y, tree[curNodeIndex].mu, tau)) {
                searchLeftSubtree(res, curNodeIndex, y, tau, displacement, d, k);
            }
            if (checkOutside(tree[curNodeIndex].vp, y, tree[curNodeIndex].mu, tau)) {
                searchRightSubtree(res, curNodeIndex, y, tau, displacement, d, k);
            }
        }
    }

    void searchRightSubtree(knnresult res, int curNodeIndex, double y, double tau, int displacement, int d, int k) {

        std::cout << "Node " << curNodeIndex
                  << " outside condition = " << checkOutside(tree[curNodeIndex].vp, y, tree[curNodeIndex].mu, tau)
                  << " -> search right subtree\n";
        curNodeIndex = moveRight(curNodeIndex);
        if (isLeaf(curNodeIndex)) {
            std::cout << "Reached node " << curNodeIndex << std::endl << std::endl;
            kNN(res, tree[curNodeIndex].data, &y, displacement, tree[curNodeIndex].len, 1, d, k);
        } else {
            if (checkInside(tree[curNodeIndex].vp, y, tree[curNodeIndex].mu, tau)) {
                searchLeftSubtree(res, curNodeIndex, y, tau, displacement, d, k);
            }
            if (checkOutside(tree[curNodeIndex].vp, y, tree[curNodeIndex].mu, tau)) {
                searchRightSubtree(res, curNodeIndex, y, tau, displacement, d, k);
            }
        }
    }

    bool checkInside(double vp, double y, double mu, double tau) {

        double dist;
        util::computeEuclideanDistance(&vp, &y, &dist, 1, 1, 1);
        if (dist < mu + tau)
            return true;

        return false;
    }

    bool checkOutside(double vp, double y, double mu, double tau) {

        double dist;
        util::computeEuclideanDistance(&vp, &y, &dist, 1, 1, 1);
        if (dist > mu - tau)
            return true;

        return false;
    }

    // int searchLeaf(double point) {                      // return the node index of the leaf that contains the point

    //     int curNodeIndex = 0;
    //     int _leftIdx;
    //     int _rightIdx;
    //     while ( !isLeaf(curNodeIndex) ) {
    //         _leftIdx = tree[curNodeIndex].leftIndex;
    //         _rightIdx = tree[curNodeIndex].rightIndex;
    //         if ( belongsToNode(point, tree[_leftIdx]) )
    //             curNodeIndex = _leftIdx;
    //         else
    //             curNodeIndex = _rightIdx;
    //     }
    //     return tree[curNodeIndex].index;
    // }

    // bool belongsToNode(double point, Node curNode) {

    //     for (int i = 0; i < curNode.len; i++) {
    //         if (curNode.data[i] == point)
    //             return true;
    //     }
    //     return false;
    // }

    // bool isLeftChild(int curNodeIndex) {

    //     if (curNodeIndex == 0)
    //         return false;
    //     if (tree[tree[curNodeIndex].parentIndex].leftIndex == curNodeIndex)
    //         return true;

    //     return false;
    // }

    // bool isRightChild(int curNodeIndex) {

    //     if (curNodeIndex == 0)
    //         return false;
    //     if (tree[tree[curNodeIndex].parentIndex].rightIndex == curNodeIndex)
    //         return true;

    //     return false;
    // }
}; // END OF VPT CLASS

double* computeLeftChild(Node* parent) {

    int _len      = (parent->len + 1) / 2;
    double* _data = new double[_len];
    int idx       = 0;
    for (int i = 0; i < parent->len; i++) {
        if (parent->D[i] <= parent->mu)
            _data[idx++] = parent->data[i];
    }
    // std::cout << "Left child:\t";
    // prt::rowMajor(_data, 1, _len);
    return _data;
}

double* computeRightChild(Node* parent) {

    int _len      = parent->len - ((parent->len + 1) / 2);
    double* _data = new double[_len];
    int idx       = 0;
    for (int i = 0; i < parent->len; i++) {
        if (parent->D[i] > parent->mu)
            _data[idx++] = parent->data[i];
    }
    // std::cout << "Right child:\t";
    // prt::rowMajor(_data, 1, _len);
    return _data;
}