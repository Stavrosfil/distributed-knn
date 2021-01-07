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

    std::cout << std::endl;

    util::Timer timer(true);

    // int d = 10;
    // int k = 100;
    // int b = 100;
    // int n = 10000;

    // std::string line;
    // std::string fileName = "data.csv";
    // std::ifstream myfile(fileName);

    // int d = 1;
    // int k = 1;
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

    int d = 10;
    int k = 50;
    int b = 50;
    int n = 90000;

    std::string line;
    std::string fileName = "data4.csv";
    std::ifstream myfile(fileName);

    if (argc == 2)
        fileName = argv[1];

    std::cout << fileName << std::endl;

    /* ----------------------------------- v1 ----------------------------------- */

    // std::vector<double> X;

    // timer.start("v1");

    // struct knnresult result = mpi::distrAllkNN(X, n, d, k, fileName);

    // timer.stop();

    // // prt::kNN(result);

    /* ----------------------------------- v2 ----------------------------------- */

    timer.start("v2");

    std::vector<double> X;

    struct knnresult result = mpi::distrVPTkNN(X, n, d, k, b, fileName);

    timer.stop();

    // prt::kNN(result);

    /* -------------------------------------------------------------------------- */

    std::cout << std::endl;

    return 0;
}