#include <iostream>
#include <stdio.h>
#include <cmath>
#include "knn.hpp"

class Node {

    public:

    int index;
    int len;
    double *data;
    double *vp;
    double *D;
    double mu;
    int parentIndex;
    int leftIndex;
    int rightIndex;


    Node() {

        index = -1;
        len = 0;
        data = NULL;
        vp = NULL;
        D = NULL;
        mu = -1;
        parentIndex = -1;
        leftIndex = -1;
        rightIndex = -1;
    }

    Node(int index,
    int len,
    double *data,
    int parentIndex) {

        this->index = index;
        this->len = len;
        this->data = data;
        this->vp = data;
        D = new double[len];
        util::computeEuclideanDistance(data, vp, D, len, 1, 1);
        this->mu = computeMu();
        this->parentIndex = parentIndex;
        this->leftIndex = -1;
        this->rightIndex = -1;
    }

    double computeMu() {   // TODO FOR MORE DIMENSIONS
    
        double *sortedD = new double[len];
        memcpy(sortedD, D, len * sizeof(double));
        qsort(sortedD, len, sizeof(double), util::compare);             // TODO more efficient
        mu = sortedD[(len - 1) / 2];
        //prt::rowMajor(sortedD, 1, len);
        delete[] sortedD;
        return mu;
    }
};  // END OF NODE CLASS

double* computeLeftChild(Node *parent) {
    
    int _len = (parent->len + 1) / 2;
    double *_data = new double[_len];
    int idx = 0;
    for (int i = 0; i < parent->len; i++) {
        if (parent->D[i] <= parent->mu)
            _data[idx++] = parent->data[i];
    }
    //std::cout << "Left child:\t"; 
    //prt::rowMajor(_data, 1, _len);
    return _data;
}

double* computeRightChild(Node *parent) {
    
    int _len = parent->len - ((parent->len + 1) / 2);
    double *_data = new double[_len];
    int idx = 0;
    for (int i = 0; i < parent->len; i++) {
        if (parent->D[i] > parent->mu)
            _data[idx++] = parent->data[i];
    }
    //std::cout << "Right child:\t"; 
    //prt::rowMajor(_data, 1, _len);  
    return _data;
}


class VPT {

    public:

    Node root;
    std::vector<Node> tree;
    
    VPT() {
        //std::cout << "VPT constructed\n";
    }

    VPT(Node root) {
        this->root = root;
        tree.push_back(this->root);
    }

    void createVPT() {

        int len = -1;
        int maxIdx = pow(2, (log2(root.len) + 1)) - 1;
        // std::cout << "Max Idx = " << maxIdx << std::endl; 
        for (int i = 0; i < maxIdx; i++) {
            len = tree[i].len;
            if (len > 1) {
                tree[i].leftIndex = tree.size();
                tree[i].rightIndex = tree.size() + 1;
                createChildren(tree[i], tree[i].len);
            }
        }
    }

    void createChildren(Node curNode, int len) {

        int leftLen = (curNode.len + 1) / 2;
        int rightLen = curNode.len - leftLen;
        Node _left = Node(tree.size(), leftLen, computeLeftChild(&curNode), curNode.index);
        tree.push_back(_left);
        Node _right = Node(tree.size(), rightLen, computeRightChild(&curNode), curNode.index);
        tree.push_back(_right);
    }

    int searchLeaf(double point) {                      // return the node index of the leaf that contains the point

        int curNodeIdx = 0;
        int _leftIdx;
        int _rightIdx;
        while ( !isLeaf(curNodeIdx) ) {
            _leftIdx = tree[curNodeIdx].leftIndex;
            _rightIdx = tree[curNodeIdx].rightIndex; 
            if ( belongsToNode(point, tree[_leftIdx]) )
                curNodeIdx = _leftIdx;
            else
                curNodeIdx = _rightIdx; 
        }
        return tree[curNodeIdx].index;
    }

    bool belongsToNode(double point, Node curNode) {

        for (int i = 0; i < curNode.len; i++) {
            if (curNode.data[i] == point)
                return true;
        }
        return false;
    }

    bool isLeaf(int nodeIdx) {
        if (tree[nodeIdx].leftIndex == -1)
            return true;
        else 
            return false;
    }

    int moveUp(int curIndex) {
        return tree[curIndex].parentIndex;
    }

    int moveLeft(int curIndex) {
        return tree[curIndex].leftIndex;
    }

    int moveRight(int curIndex) {
        return tree[curIndex].rightIndex;
    }

    // void vptKnn(knnresult res, double* X, double* y, int displacement, int n, int m, int d, int k) {

    //     int leafIndex = searchLeaf(*y);
    //     kNN(res, tree[leafIndex].data, y, displacement, tree[leafIndex].len, 1, 1, 1);
    // }
};  // END OF VPT CLASS