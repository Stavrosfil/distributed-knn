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
            node.leafPoints    = new Point[n];
            node.leafPointsLen = n;
            for (int i = 0; i < n; i++) {
                node.leafPoints[i] = _points[lo + 1 + i];
            }
        }

        _nodes.push_back(node);
        return _nodes.size() - 1;
    }

    void kNN(Point& p, knnresult& ans)
    {
        k             = ans.k;
        int leafIndex = searchLeaf(p, _nodes.size() - 1);
        leafKNN(p, _nodes[leafIndex]);
        climbVPT(p, leafIndex);

        for (int i = 0; i < k; i++) {
            ans.ndist[k - i - 1] = _heap.top().first;
            ans.nidx[k - i - 1]  = _heap.top().second.index;
            _heap.pop();
        }
    }

  private:
    int k = 3;
    int b = 0;

    pointHeap _heap;
    std::vector<Point>& _points;

    int searchLeaf(Point& p, int rootIndex)
    {
        int curNodeIndex = rootIndex;
        while (!isLeaf(_nodes[curNodeIndex])) {

            updateKNN(_heap, p, _points[_nodes[curNodeIndex].vpIndex], k);

            double dist = util::distance(p, _points[_nodes[curNodeIndex].vpIndex]);
            if (dist < _nodes[curNodeIndex].mu)
                curNodeIndex = moveLeft(curNodeIndex);
            else
                curNodeIndex = moveRight(curNodeIndex);
        }

        updateKNN(_heap, p, _points[_nodes[curNodeIndex].vpIndex], k);
        return curNodeIndex;
    }

    double computeMu(int lo, int median, int hi)
    {
        // prt::leafPoints(_points);
        std::nth_element(_points.begin() + lo + 1,
                         _points.begin() + median,
                         _points.begin() + hi,
                         comp::distanceFromVP(_points[lo]));
        // prt::leafPoints(_points);
        return util::distance(_points[lo], _points[median]);
    }

    void climbVPT(Point& p, int curNodeIndex)
    {
        do {
            if (isLeftChild(curNodeIndex)) {
                curNodeIndex = moveUp(curNodeIndex);
                if (checkOutside(p, _nodes[curNodeIndex], _heap.top().first)) {
                    searchSubtree(p, moveRight(curNodeIndex));
                }
            }
            else if (isRightChild(curNodeIndex)) {
                curNodeIndex = moveUp(curNodeIndex);
                if (checkInside(p, _nodes[curNodeIndex], _heap.top().first)) {
                    searchSubtree(p, moveLeft(curNodeIndex));
                }
            }
        } while (curNodeIndex != _nodes.size() - 1);
    }

    void leafKNN(Point& p, Node& curNode)
    {
        for (int i = 0; i < curNode.leafPointsLen; i++) {
            updateKNN(_heap, p, curNode.leafPoints[i], k);
        }
    }

    void searchSubtree(Point& p, int curNodeIndex)
    {
        updateKNN(_heap, p, _points[_nodes[curNodeIndex].vpIndex], k);

        if (isLeaf(_nodes[curNodeIndex]))
            leafKNN(p, _nodes[curNodeIndex]);
        else {
            if (checkInside(p, _nodes[curNodeIndex], _heap.top().first)) {
                searchSubtree(p, moveLeft(curNodeIndex));
            }
            if (checkOutside(p, _nodes[curNodeIndex], _heap.top().first)) {
                searchSubtree(p, moveRight(curNodeIndex));
            }
        }
    }

    /* -------------------------------- Checkers -------------------------------- */

    bool isLeaf(Node& n)
    {
        return n.leftIndex == -1 && n.rightIndex == -1;
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
        // std::cout << "vp = " << _points[curNode.vpIndex].coords[0]
        //   << "\tcheckInside = 1\t\tdistance = " << util::distance(p, _points[curNode.vpIndex])
        //   << "\tmu = " << curNode.mu << "\ttau = " << tau << std::endl;
        return util::distance(p, _points[curNode.vpIndex]) < curNode.mu + tau;
    }

    bool checkOutside(Point& p, Node& curNode, double tau)
    {
        // std::cout << "vp = " << _points[curNode.vpIndex].coords[0]
        //   << "\tcheckOutside = 0\tdistance = " << util::distance(p, _points[curNode.vpIndex])
        //   << "\tmu = " << curNode.mu << "\ttau = " << tau << std::endl;
        return util::distance(p, _points[curNode.vpIndex]) > curNode.mu - tau;
    }

    bool isLeftChild(int curNodeIndex)
    {
        return _nodes[curNodeIndex].parentIndex != -1
                   ? curNodeIndex == _nodes[_nodes[curNodeIndex].parentIndex].leftIndex
                   : false;
    }

    bool isRightChild(int curNodeIndex)
    {
        return _nodes[curNodeIndex].parentIndex != -1
                   ? curNodeIndex == _nodes[_nodes[curNodeIndex].parentIndex].rightIndex
                   : false;
    }
};

#endif // __VPT_H__