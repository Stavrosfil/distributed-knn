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
#include "distributedVPT.hpp"

int main(int argc, char** argv)
{

    /* ---------------------------- File preparations --------------------------- */

    std::cout << std::endl;

    // int d = 10;
    // int k = 100;
    // int b = 100;
    // int n = 10000;

    // std::string line;
    // std::string fileName = "data.csv";
    // std::ifstream myfile(fileName);

    // int d = 1;
    // int k = 2;
    // int b = 0;
    // int n = 8;

    // std::string line;
    // std::string fileName = "data2.csv";
    // std::ifstream myfile(fileName);

    // int d = 2;
    // int k = 2;
    // int b = 0;
    // int n = 4;

    // std::string line;
    // std::string fileName = "data3.csv";
    // std::ifstream myfile(fileName);

    int d = 10;
    int k = 50;
    int b = 50;
    int n = 30000;

    std::string line;
    std::string fileName = "data4.csv";
    std::ifstream myfile(fileName);


    if (argc == 2)
        fileName = argv[1];

    std::cout << fileName << std::endl;

    /* ----------------------------------- v1 ----------------------------------- */

    // std::vector<double> X;

    // cl::start();

    // struct knnresult result = mpi::distrAllkNN(X, n, d, k, fileName);

    // cl::stop(true, "v1");

    // // prt::kNN(result);

    /* -------------------------------- Build VPT ------------------------------- */

    // std::vector<Point> corpus(n);
    // util::readNDimVector(corpus, n, d, fileName);
    // std::vector<Point> query(corpus);

    // VPT vpt(corpus, b, k);

    // cl::start();

    // Node* root = vpt.buildTree(0, corpus.size());

    // cl::stop(true, "Vantage point tree building");

    // // std::cout << "\nCorpus points:\n";
    // // prt::points(corpus);

    // // std::cout << "\nVantage point tree:\n\n";
    // // prt::tree(root, corpus);
    // // prt::point(root->right->right->leafPoints[0]);
    // // prt::point(root->right->leafPoints[0]);

    // // std::cout << std::endl;

    /* ----------------------------------- v2 ----------------------------------- */

    // knnresult ans = knnresult();
    // ans.m         = query.size();
    // ans.k         = k;
    // ans.nidx      = new int[ans.m * ans.k];
    // ans.ndist     = new double[ans.m * ans.k];

    // std::fill_n(ans.nidx, ans.m * ans.k, D_MAX);
    // std::fill_n(ans.ndist, ans.m * ans.k, -1);

    // cl::start();

    // for (auto p : query) {
    //     vpt.kNN(p, ans, p.index, *root);
    // }

    // cl::stop(true, "v2");

    // // prt::kNN(ans);

    /* ----------------------------- Reconstruct VPT ---------------------------- */

    // VPT vpt2(corpus, b, k);

    // cl::start();

    // Node* root2 = vpt2.reconstructTree(0, corpus.size());

    // cl::stop(true, "Vantage point tree reconstruction");

    // // std::cout << "\nReconstructed vantage point tree:\n\n";

    // // prt::tree(root2, corpus);

    /* ---------------------------- Serialize Vector ---------------------------- */

    // cl::start();

    // int* indices   = new int[corpus.size()];
    // double* coords = new double[corpus[0].d * corpus.size()];

    // conv::serVector(corpus, indices, coords);

    // cl::stop(true, "Vector serialization");

    // // std::cout << "\nVector serialization: \n" << "Coords:\t\t";
    // // prt::rowMajor(coords, 1, corpus[0].d * corpus.size());
    // // std::cout << "Indices:\t";
    // // prt::rowMajor(indices, 1, corpus.size());

    /* --------------------------- Reconstruct Vector --------------------------- */

    // cl::start();

    // std::vector<Point> rV(n);
    // conv::recVector(rV, indices, coords, d);

    // cl::stop(true, "Vector reconstruction");

    // // std::cout << "\nVector reconstruction: \n";
    // // prt::points(rV);

    /* -------------------------- Distributed VPT kNN --------------------------- */

    cl::start();

    std::vector<double> X;

    struct knnresult result = mpi::distrVPTkNN(X, n, d, k, b, fileName);

    cl::stop(true, "Distributed VPT kNN");

    // prt::kNN(result);

    /* -------------------------------------------------------------------------- */

    std::cout << std::endl;

    return 0;
}