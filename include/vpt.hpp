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
    VPT(std::vector<Point>& points, int b, int k, int queryIndexNum) : _points(points), _b(b), _k(k)
    {
        _heaps.resize(queryIndexNum);
    }

    void computeKNN(Point& qp, knnresult& ans, int queryIndex)
    {
        _nodes_visits = 0;

        _heap = _heaps[queryIndex];

        _tau = D_MAX;

        search(qp, *_root);

        _heaps[queryIndex] = _heap;

        for (int i = 0; i < _k; i++) {
            int index = _k * (queryIndex + 1) - i - 1;

            ans.ndist[index] = _heap.top().dist;
            ans.nidx[index]  = _heap.top().index;

            _heap.pop();
        }
    }

    void _computeKNN(Point& qp, knnresult& ans, int queryIndex) 
    {
        _nodes_visits = 0;

        _heap = _heaps[queryIndex];

        _tau = D_MAX;

        _search(qp, *_root);

        _heaps[queryIndex] = _heap;

        for (int i = 0; i < _k; i++) {
            int index = _k * (queryIndex + 1) - i - 1;

            ans.ndist[index] = _heap.top().dist;
            ans.nidx[index]  = _heap.top().index;

            _heap.pop();
        }
    }

    void build(int points_size)
    {
        _root = buildTree(0, points_size);
    }

    void rebuild(int points_size)
    {
        _root = rebuildTree(0, points_size);
    }

    int getNodesVisits()
    {
        return _nodes_visits;
    }

  private:
    // Dynamic size
    std::vector<Point>& _points;

    // Constant size
    Node* _root;
    int _k;
    int _b;
    int _nodes_visits;
    double _tau = D_MAX;

    pointHeap _heap;
    std::vector<pointHeap> _heaps;

    // Partition and construct tree from _points[lo:hi]
    // lo and hi are indices of _points vector
    Node* buildTree(int lo, int hi)
    {
        Node* node    = new Node();
        node->vpIndex = lo;

        int n = hi - (lo + 1);
        if (n >= 2 * _b + 2) {
            int median  = (hi + lo) % 2 == 0 ? (hi + lo) / 2 : (hi + lo + 1) / 2;
            node->mu    = computeMu(lo, median, hi);
            node->left  = buildTree(lo + 1, median);
            node->right = buildTree(median, hi);
        }
        else {
            // int i = (int)((double)rand() / RAND_MAX * (hi - lo - 1)) + lo;
            // std::swap(_points[lo], _points[i]);

            node->leafPointsIndex = lo + 1;
            node->leafPointsLen   = n;
        }

        return node;
    }    

    Node* rebuildTree(int lo, int hi)
    {
        Node* node    = new Node();
        node->vpIndex = lo;

        int n = hi - (lo + 1);
        if (n >= 2 * _b + 2) {
            int median  = (hi + lo) % 2 == 0 ? (hi + lo) / 2 : (hi + lo + 1) / 2;
            node->mu    = util::distance(_points[lo], _points[median]);
            node->left  = rebuildTree(lo + 1, median);
            node->right = rebuildTree(median, hi);
        }
        else {
            // int i = (int)((double)rand() / RAND_MAX * (hi - lo - 1)) + lo;
            // std::swap(_points[lo], _points[i]);

            node->leafPointsIndex = lo + 1;
            node->leafPointsLen   = n;
        }

        return node;
    }

    void search(Point& qp, Node& node)
    {
        if (node.vpIndex == -1)
            return;

        _nodes_visits++;

        double dist = util::distance(qp, _points[node.vpIndex]);

        if (dist < _tau) {
            if (_heap.size() == _k)
                _heap.pop();
            _heap.push(heapItem(dist, _points[node.vpIndex].index));
            if (_heap.size() == _k)
                _tau = _heap.top().dist;
        }

        if (isLeaf(node)) {
            leafKNN(qp, node);
            if (_heap.size() == _k)
                _tau = _heap.top().dist;
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

    void _search(Point& qp, Node& root) 
    {
        // if (root.vpIndex == -1)
        //     return;

        _nodes_visits++;

        // kNN with vp
        double dist = updateKNN(_heap, qp, _points[root.vpIndex], _k);
        
        // update tau
        if (_heap.size() == _k)
            _tau = _heap.top().dist;

        // leaf kNN & update tau
        if (isLeaf(root)) {
            leafKNN(qp, root);
            if (_heap.size() == _k)
                _tau = _heap.top().dist;
            return;
        }

        // search
        if (dist < root.mu) {
            _search(qp, *root.left);
            
            if (_checkRightInt(qp, root)) {
                _searchSubtree(qp, *root.right);
            }
        }
        else {
            _search(qp, *root.right);

            if (_checkLeftInt(qp, root)) {
                _searchSubtree(qp, *root.left);
            }
        }
    }

    void _searchSubtree(Point& qp, Node& subroot) 
    {
        // kNN with vp
        double dist = updateKNN(_heap, qp, _points[subroot.vpIndex], _k);
        
        // update tau
        if (_heap.size() == _k)
            _tau = _heap.top().dist;

        // leaf kNN & update tau
        if (isLeaf(subroot)) {
            leafKNN(qp, subroot);
            if (_heap.size() == _k)
                _tau = _heap.top().dist;
            return;
        }

        // search
        if (_checkLeftInt(qp, subroot)) {
            _searchSubtree(qp, *subroot.left);
        }

        if (_checkRightInt(qp, subroot)) {
            _searchSubtree(qp, *subroot.right);
        }
    }

    // TODO optimize leafKNN, use euclideanDistance vector
    void leafKNN(Point& qp, Node& curNode)
    {
        for (int i = 0; i < curNode.leafPointsLen; i++) {
            _nodes_visits++;
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

    bool _checkLeftInt(Point& qp, Node& node)
    {
        return util::distance(qp, _points[node.vpIndex]) < node.mu + _tau;
    }

    bool _checkRightInt(Point& qp, Node& node)
    {
        return util::distance(qp, _points[node.vpIndex]) > node.mu - _tau;
    }
};

#endif // __VPT_H__