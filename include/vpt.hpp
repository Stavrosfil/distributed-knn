#ifndef __VPT_H__
#define __VPT_H__

#include <array>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <queue>

#include "knn.hpp"

class VPT {

public:
    VPT(std::vector<Point>& points, int b, int k) : _points(points), _b(b), _k(k) {}

    // Partition and construct tree from _points[lo:hi]
    // lo and hi are indices of corpus array
    Node* buildTree(int lo, int hi)
    {
        Node* node = new Node();
        node->vpIndex = lo;

        int n = hi - (lo + 1);
        if (n >= 2 * _b + 2) {
            int median = (hi + lo) % 2 == 0 ? (hi + lo) / 2 : (hi + lo + 1) / 2;
            node->mu = computeMu(lo, median, hi);
            node->left = buildTree(lo + 1, median);
            node->right = buildTree(median, hi);
        }
        else {
            // int i = (int)((double)rand() / RAND_MAX * (hi - lo - 1)) + lo;
            // std::swap(_points[lo], _points[i]);
            // std::cout << "leafPointsLen = " << n << std::endl;
            node->leafPointsIndex = lo + 1;
            node->leafPointsLen = n;
 
            // for (int i = 0; i < node->leafPointsLen; i++) {
            //     prt::point(_points[node->leafPointsIndex + i]);
            // }
        }

        return node;
    }

    void kNN(Point& qp, knnresult& ans, int queryIndex, Node& root)
    {
        _heap.push(heapItem(D_MAX, nullptr));

        _tau = D_MAX;

        search(qp, root);
       
        for (int i = 0; i < _k; i++) {
            ans.ndist[_k * (queryIndex + 1) - i - 1] = _heap.top().dist;
            ans.nidx[_k * (queryIndex + 1) - i - 1] = _heap.top().p->index;
            _heap.pop();
        }
    }

    Node* reconstructTree(int lo, int hi) 
    {
        Node* node = new Node();
        node->vpIndex = lo;

        int n = hi - (lo + 1);
        if (n >= 2 * _b + 2) {
            int median = (hi + lo) % 2 == 0 ? (hi + lo) / 2 : (hi + lo + 1) / 2;
            node->mu = util::distance(_points[lo], _points[median]);
            node->left = reconstructTree(lo + 1, median);
            node->right = reconstructTree(median, hi);
        }
        else {
            // int i = (int)((double)rand() / RAND_MAX * (hi - lo - 1)) + lo;
            // std::swap(_points[lo], _points[i]);

            node->leafPointsIndex = lo + 1;
            node->leafPointsLen = n;
 
            // for (int i = 0; i < node->leafPointsLen; i++) {
            //     prt::point(_points[node->leafPointsIndex + i]);
            // }
        }

        return node;
    }

private:
    int _k;
    int _b;
    double _tau = D_MAX;
    std::vector<Point>& _points;

    pointHeap _heap;

    void search(Point& qp, Node& node)
    {
        if (node.vpIndex == -1)
            return;

        double dist = util::distance(qp, _points[node.vpIndex]);

        if (dist < _tau) {
            if (_heap.size() == _k)
                _heap.pop();
            _heap.push(heapItem(dist, &_points[node.vpIndex]));
            if (_heap.size() == _k)
                _tau = _heap.top().dist;
        }

        if (isLeaf(node)) {
            leafKNN(qp, node);
            return;
        }

        if (dist < node.mu) {
            if (dist - _tau <= node.mu)
                search(qp, *node.left);
            if (dist + _tau >= node.mu)
                search(qp, *node.right);
        }
        else {
            if (dist + _tau >= node.mu)
                search(qp, *node.right);
            if (dist - _tau <= node.mu)
                search(qp, *node.left);
        }
    }

    // TODO optimize leafKNN
    void leafKNN(Point& qp, Node& curNode)
    {
        for (int i = 0; i < curNode.leafPointsLen; i++) {
            // std::cout << "LeafKNN\n" << curNode.leafPointsLen;
            updateKNN(_heap, qp, _points[curNode.leafPointsIndex + i], _k);
        }
    }

    double computeMu(int lo, int median, int hi)
    {
        std::nth_element(_points.begin() + lo + 1,
            _points.begin() + median,
            _points.begin() + hi,
            comp::distanceFromVP(_points[lo]));
        return util::distance(_points[lo], _points[median]);
    }

    bool isLeaf(Node& n)
    {
        return !n.left && !n.right;
    }
};

#endif // __VPT_H__