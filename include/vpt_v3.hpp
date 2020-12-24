#include <iostream>
#include <stdio.h>
#include <cmath>

class Node {

    public:

    int index;
    int len;
    double *data;
    double *vp;
    double *D;
    double mu;
    // Node *parent;
    // Node *left;
    // Node *right;
    int parentIndex;
    int leftIndex;
    int rightIndex;


    Node() {

        //std::cout << "Node constructed\n";
        index = -1;
        len = 0;
        data = NULL;
        vp = NULL;
        D = NULL;
        mu = -1;
        // parent = NULL;
        // left = NULL;
        // right = NULL;
        parentIndex = -1;
        leftIndex = -1;
        rightIndex = -1;
    }

    Node(int index,
    int len,
    double *data,
    int parentIndex) {

        //std::cout << "Node constructed\n";
        this->index = index;
        this->len = len;
        this->data = data;
        this->vp = data;
        D = new double[len];
        util::computeEuclideanDistance(data, vp, D, len, 1, 1);
        this->mu = computeMu();
        this->parentIndex = parentIndex;
        this->leftIndex = index + 1;
        this->rightIndex = index + 2;
    }

    double computeMu() {   // TODO FOR MORE DIMENSIONS
    
        double *sortedD = new double[len];
        memcpy(sortedD, D, len * sizeof(double));
        qsort(sortedD, len, sizeof(double), util::compare);             // TODO more efficient
        mu = sortedD[len / 2];
        //prt::rowMajor(sortedD, 1, len);
        return mu;
    }
};

// void computeLeft(Node *node) {

//         double *_data = new double[node->len / 2 + 1];
//         int idx = 0;
//         if (node->len != 0) {
//             for (int i = 0; i < node->len; i++) {
//                 if (node->D[i] <= node->mu)
//                     _data[idx++] = node->data[i];
//             }

//             //Node temp = Node(node->index + 1, node->len / 2 + 1, _data, NULL);
//             node->left = Node(node->index + 1, node->len / 2 + 1, _data, NULL);
//         }
// }

double* computeLeftChild(Node *parent) {
    
    int _len = parent->len / 2 + 1;
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
    
    int _len = parent->len - (parent->len / 2 + 1);
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
        std::cout << "VPT constructed\n";
    }

    VPT(Node root) {
        this->root = root;
        tree.push_back(this->root);
    }

    void createVPT() {

        if (root.len > 1) {
            int leftLen = root.len / 2 + 1;
            int rightLen = root.len - leftLen;
            Node left = Node(root.index + 1,leftLen, computeLeftChild(&root), -1);
            tree.push_back(left);
            Node right = Node(root.index + 2, rightLen, computeRightChild(&root), -1);
            tree.push_back(right);
        }

    }
};