#include <iostream>
#include <stdlib.h>
#include <iterator>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <chrono>
// #include <benchmark/benchmark.h>

#include "knn.hpp"
#include "utils.hpp"
#include "distributed.hpp"
#include "vpt.hpp"
#include "distributedVPT.hpp"

int main(int argc, char** argv)
{

    /* ---------------------------- File preparations --------------------------- */

    // int d = 1;
    // int k = 8;
    // int b = 0;
    // int n = 8;

    // std::string fileName = "data2.csv";

    // int d = 2;
    // int k = 4;
    // int b = 0;
    // int n = 4;

    // std::string fileName = "data3.csv";

    int d = 32;
    int k = 50;
    int b = 50;
    int n = 68040;

    std::string fileName = "ColorHistogram.asc";

    if (argc == 2)
        fileName = argv[1];

    /* ----------------------------------- v1 ----------------------------------- */

    // std::vector<double> X;

    // struct knnresult result = mpi::distrAllkNN(X, n, d, k, fileName);

    // // prt::kNN(result);

    /* ----------------------------------- v2 ----------------------------------- */

    std::vector<double> X;

    struct knnresult result = mpi::distrVPTkNN(X, n, d, k, b, fileName);

    // prt::kNN(result);

    /* -------------------------------------------------------------------------- */

    return 0;
}