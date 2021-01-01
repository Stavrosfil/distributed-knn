#ifndef __VPT_H__
#define __VPT_H__

#include <array>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <queue>

#include "knn.hpp"

struct distanceFromVP {
    Point& p;
    distanceFromVP(Point& p) : p(p) {}
    bool operator()(Point& p1, Point& p2)
    {
        return util::distance(p, p1) < util::distance(p, p2);
    }
};

class VPT {

  public:
    std::vector<Node> _nodes;

    VPT(std::vector<Point>& points) : _points(points) {}

    // Partition and construct tree from _points[lo:hi]
    // lo and hi are indices of corpus array
    int buildTree(int lo, int hi)
    {
        Node node;
        node.vpIndex = lo;

        int n = hi - (lo + 1);
        if (n >= 2 * b + 2) {
            int median                          = (hi + lo) % 2 == 0 ? (hi + lo) / 2 : (hi + lo + 1) / 2;
            node.mu                             = computeMu(lo, median, hi);
            node.leftIndex                      = buildTree(lo + 1, median);
            node.rightIndex                     = buildTree(median, hi);
            _nodes[node.leftIndex].parentIndex  = _nodes.size();
            _nodes[node.rightIndex].parentIndex = _nodes.size();
        }
        else {
            node.points     = new Point[n];
            node.points_len = n;
            for (int i = 0; i < n; i++) {
                node.points[i] = _points[lo + 1 + i];
            }
        }

        _nodes.push_back(node);
        return _nodes.size() - 1;
    }

    void vptKnn(Point& p, knnresult& ans)
    {
        // heap.push(HeapItem(_points[0], util::distance(p, _points[0])));
        k             = ans.k;
        int leafIndex = searchLeaf(p, _nodes.size() - 1);
        leafKNN(p, _nodes[leafIndex]);
        climbVPT(p, leafIndex);

        for (int i = 0; i < k; i++) {
            ans.ndist[k - i - 1] = heap.top().first;
            ans.nidx[k - i - 1]  = heap.top().second.index;
            heap.pop();
        }
    }

  private:
    int k = 3;
    int b = 0;
    double tau;
    knnresult& ans();
    std::vector<Point>& _points;

    std::priority_queue<std::pair<double, Point>, std::vector<std::pair<double, Point>>, CustomCompare> heap;

    int searchLeaf(Point& p, int rootIndex)
    {
        int curNodeIndex = rootIndex;
        while (!isLeaf(_nodes[curNodeIndex])) {

            updateKNN(heap, p, _points[_nodes[curNodeIndex].vpIndex], k);

            double dist = util::distance(p, _points[_nodes[curNodeIndex].vpIndex]);
            if (dist < _nodes[curNodeIndex].mu)
                curNodeIndex = moveLeft(curNodeIndex);
            else
                curNodeIndex = moveRight(curNodeIndex);
        }

        updateKNN(heap, p, _points[_nodes[curNodeIndex].vpIndex], k);

        std::cout << "Leaf vp = " << _points[_nodes[curNodeIndex].vpIndex].coords[0] << std::endl;
        return curNodeIndex;
    }

    double computeMu(int lo, int median, int hi)
    {
        // prt::points(_points);
        std::nth_element(
            _points.begin() + lo + 1, _points.begin() + median, _points.begin() + hi, distanceFromVP(_points[lo]));
        // prt::points(_points);
        return util::distance(_points[lo], _points[median]);
    }

    void climbVPT(Point& p, int curNodeIndex)
    { // curNodeIndex = leaf index
        do {
            if (isLeftChild(curNodeIndex)) {
                curNodeIndex = moveUp(curNodeIndex);
                if (checkOutside(p, _nodes[curNodeIndex], heap.top().first)) {
                    searchSubtree(p, moveRight(curNodeIndex));
                }
            }
            else if (isRightChild(curNodeIndex)) {
                curNodeIndex = moveUp(curNodeIndex);
                if (checkInside(p, _nodes[curNodeIndex], heap.top().first)) {
                    searchSubtree(p, moveLeft(curNodeIndex));
                }
            }
        } while (curNodeIndex != _nodes.size() - 1);
    }

    void leafKNN(Point& p, Node& curNode)
    {
        if (isLeaf(curNode)) { // TODO remove if
            for (int i = 0; i < curNode.points_len; i++) {
                updateKNN(heap, p, curNode.points[i], k);
            }
        }
    }

    bool isLeaf(Node& n)
    {
        if (n.leftIndex == -1)
            return true;
        else
            return false;
    }

    int moveUp(int curNodeIndex)
    {
        return _nodes[curNodeIndex].parentIndex;
    }

    int moveLeft(int curNodeIndex)
    {
        return _nodes[curNodeIndex].leftIndex;
    }

    int moveRight(int curNodeIndex)
    {
        return _nodes[curNodeIndex].rightIndex;
    }

    bool checkInside(Point& p, Node& curNode, double tau)
    {
        if (util::distance(p, _points[curNode.vpIndex]) < curNode.mu + tau) {
            std::cout << "vp = " << _points[curNode.vpIndex].coords[0]
                      << "\tcheckInside = 1\t\tdistance = " << util::distance(p, _points[curNode.vpIndex])
                      << "\tmu = " << curNode.mu << "\ttau = " << tau << std::endl;
            return true;
        }
        std::cout << "vp = " << _points[curNode.vpIndex].coords[0]
                  << "\tcheckInside = 0\t\tdistance = " << util::distance(p, _points[curNode.vpIndex])
                  << "\tmu = " << curNode.mu << "\ttau = " << tau << std::endl;
        return false;
    }

    bool checkOutside(Point& p, Node& curNode, double tau)
    {
        if (util::distance(p, _points[curNode.vpIndex]) > curNode.mu - tau) {
            std::cout << "vp = " << _points[curNode.vpIndex].coords[0]
                      << "\tcheckOutside = 1\tdistance = " << util::distance(p, _points[curNode.vpIndex])
                      << "\tmu = " << curNode.mu << "\ttau = " << tau << std::endl;
            return true;
        }
        std::cout << "vp = " << _points[curNode.vpIndex].coords[0]
                  << "\tcheckOutside = 0\tdistance = " << util::distance(p, _points[curNode.vpIndex])
                  << "\tmu = " << curNode.mu << "\ttau = " << tau << std::endl;
        return false;
    }

    bool isLeftChild(int curNodeIndex)
    {
        if (_nodes[curNodeIndex].parentIndex == -1)
            return false;
        if (curNodeIndex == _nodes[_nodes[curNodeIndex].parentIndex].leftIndex)
            return true;
        return false;
    }

    bool isRightChild(int curNodeIndex)
    {
        if (_nodes[curNodeIndex].parentIndex == -1)
            return false;
        if (curNodeIndex == _nodes[_nodes[curNodeIndex].parentIndex].rightIndex)
            return true;
        return false;
    }

    void searchSubtree(Point& p, int curNodeIndex)
    {
        updateKNN(heap, p, _points[_nodes[curNodeIndex].vpIndex], k);

        if (isLeaf(_nodes[curNodeIndex]))
            leafKNN(p, _nodes[curNodeIndex]);
        else {
            if (checkInside(p, _nodes[curNodeIndex], heap.top().first)) {
                searchSubtree(p, moveLeft(curNodeIndex));
            }
            if (checkOutside(p, _nodes[curNodeIndex], heap.top().first)) {
                searchSubtree(p, moveRight(curNodeIndex));
            }
        }
    }
};

#endif // __VPT_H__