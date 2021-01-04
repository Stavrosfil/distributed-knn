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
#include "vptree.hpp"

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

    int d = 10;
    int k = 100;
    int b = 0;
    int n = 10000;

    std::string line;
    std::string fileName = "data.csv";
    std::ifstream myfile(fileName);

    // int d = 1;
    // int k = 4;
    // int b = 0;
    // int n = 8;

    // std::string line;
    // std::string fileName = "data2.csv";
    // std::ifstream myfile(fileName);

    // int d = 2;
    // int k = 4;
    // int b = 0;
    // int n = 4;

    // std::string line;
    // std::string fileName = "data3.csv";
    // std::ifstream myfile(fileName);

    if (argc == 2)
        fileName = argv[1];

    std::cout << fileName << std::endl;

    /* ------------------------------- Init clock ------------------------------- */

    auto t1 = std::chrono::high_resolution_clock::now();
    auto t2 = std::chrono::high_resolution_clock::now();
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

    /* -------------------------------- Build VPT ------------------------------- */

    std::vector<Point> corpus = getNDVector(n, d, fileName);
    std::vector<Point> query(corpus);

    VPT vpt(corpus, b, k);             // TODO pass corpus by reference and make _points private

    t1 = std::chrono::high_resolution_clock::now();

    Node* root = vpt.buildTree(0, corpus.size());

    t2 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

    std::cout << "\nVantage point tree building time: " << duration / 1e3 << "ms" << std::endl;

    // prt::points(corpus);

    /* ----------------------------------- v2 ----------------------------------- */

    knnresult ans = knnresult();
    ans.m = query.size();
    ans.k = k;
    ans.nidx = new int[ans.m * ans.k];
    ans.ndist = new double[ans.m * ans.k];

    std::fill_n(ans.nidx, ans.m * ans.k, D_MAX);
    std::fill_n(ans.ndist, ans.m * ans.k, -1);

    // t1 = std::chrono::high_resolution_clock::now();

    // for (auto p : query) {
    //     vpt.kNN(p, ans, p.index, *root);
    // }

    // t2 = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

    // prt::kNN(ans);

    // std::cout << duration / 1e3 << "ms" << std::endl;

    /* --------------------------------- Prints --------------------------------- */

    // std::cout << "\nVantage point tree:\n\n";

    // prt::tree(root, vpt._points);

    // // prt::point(root->right->right->leafPoints[0]);

    // std::cout << std::endl;

    /* ----------------------------- Reconstruct VPT ---------------------------- */

    VPT vpt2(corpus, b, k);
    
    t1 = std::chrono::high_resolution_clock::now();

    Node* root2 = vpt2.reconstructTree(0, corpus.size());
    
    t2 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

    std::cout << "\nVantage point tree reconstruction time: " << duration / 1e3 << "ms" << std::endl;

    // std::cout << "\nReconstructed vantage point tree:\n\n";

    // prt::tree(root2, corpus);

    // prt::point(root2->right->right->leafPoints[0]);

    /* ---------------------------- Serialize Vector ---------------------------- */

    t1 = std::chrono::high_resolution_clock::now();

    int *indices = new int[corpus.size()];
    double* coords = new double[corpus[0].d * corpus.size()];

    conv::serVector(corpus, indices, coords);
    
    t2 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

    std::cout << "\nVector serialization time: " << duration / 1e3 << "ms" << std::endl;
    
    // std::cout << "\nVector serialization: \n" << "Coords:\t\t";
    // prt::rowMajor(coords, 1, corpus[0].d * corpus.size());
    // std::cout << "Indices:\t";
    // prt::rowMajor(indices, 1, corpus.size());

    /* --------------------------- Reconstruct Vector --------------------------- */

    t1 = std::chrono::high_resolution_clock::now();
    
    std::vector<Point> rV(n);
    conv::recVector(rV, indices, coords, n, d);
    
    t2 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

    std::cout << "\nVector reconstruction time: " << duration / 1e3 << "ms" << std::endl;

    // std::cout << "\nVector reconstruction: \n";
    // prt::points(rV);

    std::cout << std::endl;

    return 0;
}