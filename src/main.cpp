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

    /* ----------------------------- Init variables ----------------------------- */

    int k = 50;

    int n       = -1;
    int d       = -1;
    int b       = -1;
    int data    = -1;
    int version = -1;

    /* dataset:
        data = 0 -> corel, colorhistogram
        data = 1 -> corel, colormoments
        data = 2 -> corel, cooctexture
        data = 3 -> fma, features
        data = 4 -> miniboone
        data = 5 -> tv_news_com, BBC
        data = 6 -> tv_news_com, CNN
        data = 7 -> tv_news_com, CNNIBN
        data = 8 -> tv_news_com, NDTV
        data = 9 -> tv_news_com, TIMESNOW
    */

    if (argc == 3) {
        data    = std::stof(argv[1]);
        version = std::stof(argv[2]);
    }

    /* ------------------------------ Reader test ------------------------------- */

    // std::vector<double> X;

    // rdMiniboone::mnbPid(n, d, X, 0);

    // // prt::vector(X);

    /* ----------------------------------- v1 ----------------------------------- */

    if (version == 1) {

        std::vector<double> X;

        struct knnresult result = mpi::distrAllkNN(X, n, d, k, data);

        // prt::kNN(result);
    }

    /* ----------------------------------- v2 ----------------------------------- */

    if (version == 2) {

        std::vector<double> X;

        struct knnresult result = mpi::distrVPTkNN(X, n, d, k, b, data);

        // prt::kNN(result);
    }

    /* -------------------------------------------------------------------------- */

    return 0;
}