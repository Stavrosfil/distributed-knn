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

    // int d = 1;
    // int k = 8;
    // int b = 0;
    // int n = 8;

    // std::string line;
    // std::string fileName = "data2.csv";
    // std::string filePath = "./dataset/random/" + fileName;
    // std::ifstream myfile(filePath);
    // std::cout << fileName << std::endl;


    // int d = 2;
    // int k = 4;
    // int b = 0;
    // int n = 4;

    // std::string line;
    // std::string fileName = "data3.csv";
    // std::string filePath = "./dataset/random/" + fileName;
    // std::ifstream myfile(filePath);
    // std::cout << fileName << std::endl;

    // int d = 10;
    // int k = 50;
    // int b = 50;
    // int n = 90000;

    // std::string line;
    // std::string fileName = "data4.csv";
    // std::string filePath = "./dataset/random/" + fileName;
    // std::ifstream myfile(filePath);
    // std::cout << fileName << std::endl;


    int d = 32;
    int k = 50;
    int b = 50;
    int n = 68040;

    std::string line;
    std::string fileName = "ColorHistogram.asc";
    std::string filePath = "./dataset/corel_image_features/data/" + fileName;
    std::ifstream myfile(filePath);
    std::cout << fileName << std::endl;

    if (argc == 2)
        fileName = argv[1];

    /* ----------------------------------- v1 ----------------------------------- */

    std::vector<double> X;

timer.start("v1");

    struct knnresult result = mpi::distrAllkNN(X, n, d, k, filePath);

timer.stop();

    // prt::kNN(result);

    /* ----------------------------------- v2 ----------------------------------- */

//     std::vector<double> X;

// timer.start("v2");

//     struct knnresult result = mpi::distrVPTkNN(X, n, d, k, b, filePath);

// timer.stop();

//     // prt::kNN(result);

    /* -------------------------------------------------------------------------- */

    std::cout << std::endl;

    return 0;
}