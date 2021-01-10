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

    int n = 0;
    int d = 0;
    int k = 50;
    int b = 50;
    std::string fileName = " ";

    // if (argc == 2)
    //     fileName = argv[1];

    /* ------------------------------ Reader test ------------------------------- */



    // std::vector<double> X;

    // rdCorel::coocTex(n, d, X, 0);

    
    // // prt::vector(X);


    /* ----------------------------------- v1 ----------------------------------- */

    // std::vector<double> X;

    // struct knnresult result = mpi::distrAllkNN(X, n, d, k, fileName);

    // // prt::kNN(result);

    /* ----------------------------------- v2 ----------------------------------- */

    // std::vector<double> X;

    // struct knnresult result = mpi::distrVPTkNN(X, n, d, k, b, fileName);

    // // prt::kNN(result);

    /* -------------------------------------------------------------------------- */

    return 0;
}