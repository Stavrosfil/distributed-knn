#include <iostream>
#include <stdio.h>
#include <cmath>
#include <array>
#include <vector>

#include "knn.hpp"


class VPT {

  public:
    std::vector<Point>& _points;
    std::vector<Node> _nodes;
    int b = 0;

    VPT(std::vector<Point>& points) : _points(points) {}

    // Partition and construct tree from _points[lo:hi]
    // lo and hi are indices of corpus array
    int buildTree(int lo, int hi) {

        Node node;
        node.vpIndex = lo;

        int n = hi - (lo + 1);
        if (n >= 2 * b + 2) {
            int median      = (hi + lo) % 2 == 0 ? (hi + lo) / 2 : (hi + lo + 1) / 2;
            node.mu         = computeMu(lo, median, hi);
            node.leftIndex  = buildTree(lo + 1, median);
            node.rightIndex = buildTree(median, hi);
            _nodes[node.leftIndex].parentIndex = _nodes.size();
            _nodes[node.rightIndex].parentIndex = _nodes.size();
        } else {
            node.points     = new Point[n];
            node.points_len = n;
            for (int i = 0; i < n; i++) {
                node.points[i] = _points[lo + 1 + i];
            }
        }

        _nodes.push_back(node);
        return _nodes.size() - 1;
    }

  // private:
    double computeMu(int lo, int median, int hi) {

        // prt::points(_points);
        std::nth_element(_points.begin() + lo + 1, _points.begin() + median, _points.begin() + hi,
                         distanceFromVP(_points[lo]));
        // prt::points(_points);
        return util::distance(_points[lo], _points[median]);
    }

    struct distanceFromVP {

        Point& p;
        distanceFromVP(Point& p) : p(p) {}
        bool operator()(Point& p1, Point& p2) {
            return util::distance(p, p1) < util::distance(p, p2);
        }
    };

    int searchLeaf(Point& p, knnresult &res) {

        int curNodeIndex = _nodes.size() - 1;                   // start from root
        while ( !isLeaf(_nodes[curNodeIndex]) ) {
            kNN(res, _points[_nodes[curNodeIndex].vpIndex].coords, p.coords, _points[_nodes[curNodeIndex].vpIndex].index, 1, 1, p.d, res.k);
            double dist = util::distance(p, _points[_nodes[curNodeIndex].vpIndex]);
            if (dist < _nodes[curNodeIndex].mu)
                curNodeIndex = moveLeft(curNodeIndex);
            else
                curNodeIndex = moveRight(curNodeIndex);
        }

        kNN(res, _points[_nodes[curNodeIndex].vpIndex].coords, p.coords, _points[_nodes[curNodeIndex].vpIndex].index, 1, 1, p.d, res.k);
        return curNodeIndex;
    }

    void leafKNN(Point& p, Node& curNode, knnresult &res) {

        if ( isLeaf(curNode) ) {
            for (int i = 0; i < curNode.points_len; i++) {
                kNN(res, curNode.points[i].coords, p.coords, curNode.points[i].index, 1, 1, p.d, res.k);
            }
        }
    }

    bool isLeaf(Node& n) {
        if (n.leftIndex == -1)
            return true;
        else 
            return false;
    }

    int moveUp(int curNodeIndex) {
        return _nodes[curNodeIndex].parentIndex;
    }

    int moveLeft(int curNodeIndex) {
        return _nodes[curNodeIndex].leftIndex;
    }

    int moveRight(int curNodeIndex) {
        return _nodes[curNodeIndex].rightIndex;
    }

    bool checkInside(Point& p, Node& curNode, double tau) {

        if (util::distance(p, _points[curNode.vpIndex]) < curNode.mu + tau) {
            std::cout << "vp = " << _points[curNode.vpIndex].coords[0] << "\tcheckInside = 1\t\tdistance = " << util::distance(p, _points[curNode.vpIndex]) << "\tmu = " << curNode.mu << "\ttau = " << tau << std::endl;
            return true;
        }
        return false;
    }

    bool checkOutside(Point& p, Node& curNode, double tau) {

        if (util::distance(p, _points[curNode.vpIndex]) > curNode.mu - tau) {
            std::cout << "vp = " << _points[curNode.vpIndex].coords[0] << "\tcheckOutside = 1\tdistance = " << util::distance(p, _points[curNode.vpIndex]) << "\tmu = " << curNode.mu << "\ttau = " << tau << std::endl;
            return true;
        }
        return false;
    }

    bool isLeftChild(int curNodeIndex) {

        if (_nodes[curNodeIndex].parentIndex == -1)
            return false;
        if (curNodeIndex == _nodes[_nodes[curNodeIndex].parentIndex].leftIndex)
            return true;
        return false;
    }

    bool isRightChild(int curNodeIndex) {

        if (_nodes[curNodeIndex].parentIndex == -1)
            return false;
        if (curNodeIndex == _nodes[_nodes[curNodeIndex].parentIndex].rightIndex)
            return true;
        return false;
    }
};
