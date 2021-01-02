#include <iostream>
#include <stdlib.h>
#include <iterator>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <chrono>

#include "knn.hpp"
#include "utils.hpp"
#include "distributed.hpp"
#include "vpt.hpp"

std::vector<Point> getNDVector(int n, int d, std::string filename)
{
    std::vector<Point> result(n);

    std::ifstream input(filename);

    std::string s;

    std::getline(input, s);
    std::istringstream iss(s);
    std::cout << s << std::endl;

    for (int i = 0; i < n; i++) {
        std::getline(input, s);
        std::istringstream iss(s);

        std::string num;
        Point p(i, new double[d], d);
        int j = 0;
        while (std::getline(iss, num, ',')) {
            p.coords[j++] = std::stof(num);
        }
        result[i] = p;
    }

    return result;
}

std::vector<double> getRowMajorVector(int n, int d, std::string filename)
{
    std::vector<double> result(n * d);

    std::ifstream input(filename);

    std::string s;

    std::getline(input, s);
    std::istringstream iss(s);
    std::cout << s << std::endl;

    for (int i = 0; i < n; i++) {
        std::getline(input, s);
        std::istringstream iss(s);

        std::string num;
        int j = 0;
        while (std::getline(iss, num, ',')) {
            result[i * d + j++] = std::stof(num);
        }
    }

    return result;
}

int main(int argc, char** argv)
{

    /* ---------------------------- File preparations --------------------------- */

    std::cout << std::endl;

    int d = 5;
    int k = 50;
    int n = 10000;

    std::string line;
    std::string fileName = "data.csv";
    std::ifstream myfile(fileName);

    if (argc == 2)
        fileName = argv[1];

    std::cout << fileName << std::endl;

    // /* ------------------------------ Construct VPT ----------------------------- */

    // for (auto i : corpus) {
    //     prt::point(i);
    //     std::cout << std::endl;
    // }

    auto t1       = std::chrono::high_resolution_clock::now();
    auto t2       = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

    /* ----------------------------------- v1 ----------------------------------- */

    // std::vector<double> v1c = getRowMajorVector(n, d, fileName);
    // // prt::rowMajor(v1c.data(), n, d);

    // t1 = std::chrono::high_resolution_clock::now();

    // struct knnresult result = mpi::distrAllkNN(v1c.data(), n, d, k);

    // t2       = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

    // // prt::kNN(result);

    // std::cout << duration / 1e3 << "ms" << std::endl;

    /* ----------------------------------- v2 ----------------------------------- */

    std::vector<Point> corpus = getNDVector(n, d, fileName);
    std::vector<Point> query(corpus);

    knnresult ans = knnresult();
    ans.m         = query.size();
    ans.k         = k;
    ans.nidx      = new int[ans.m * ans.k];
    ans.ndist     = new double[ans.m * ans.k];

    std::fill_n(ans.nidx, ans.m * ans.k, D_MAX);
    std::fill_n(ans.ndist, ans.m * ans.k, -1);

    VPT vpt(corpus);
    vpt.buildTree(0, corpus.size());
    t1 = std::chrono::high_resolution_clock::now();

    // for (auto p : query)
    auto p = query[0];
    vpt.kNN(p, ans, p.index);

    t2       = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

    // prt::kNN(ans);

    std::cout << duration / 1e3 << "ms" << std::endl;

    /* --------------------------------- Prints --------------------------------- */

    // int cnt = 0;
    // std::cout << "\nVantage point tree:\n\n";
    // for (auto p : vpt._nodes) {
    //     std::cout << "nodes[" << cnt++ << "]:\tvp = ";
    //     prt::point(corpus[p.vpIndex]);
    //     std::cout << "\t\t"
    //               << "mu = " << p.mu << "\t\t(" << p.leftIndex << ", " << p.rightIndex << ", " << p.parentIndex <<
    //               ")"
    //               << std::endl;
    //     if (p.leafPointsLen) {
    //         for (int j = 0; j < p.leafPointsLen; j++) {
    //             std::cout << p.leafPoints[j].coords[0] << std::endl;
    //         }
    //     }
    // }
    // std::cout << std::endl;

    // /* --------------------------------- VPT kNN -------------------------------- */

    // for (auto p : query) {
    //     vpt.kNN(p, ans, p.index);
    // }

    // prt::kNN(ans);

    return 0;
}
