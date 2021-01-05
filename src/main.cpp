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
    // int k = 8;
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
    int n = 15000;

    std::string line;
    std::string fileName = "data4.csv";
    std::ifstream myfile(fileName);


    if (argc == 2)
        fileName = argv[1];

    std::cout << fileName << std::endl;

    /* ----------------------------------- v1 ----------------------------------- */

    // // std::vector<double> v1c(n * d);
    // // util::readToRowMajorVector(v1c, n, d, fileName);
    // // // prt::rowMajor(v1c.data(), n, d);

    // // double* X;
    // // X = new double[n * d];
    // // X = v1c.data();
    // // prt::rowMajor(X, n, d);

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

    // std::cout << std::endl;

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

    /* -------------------------------------------------------------------------- */
    



    /* ----------------------------- Distributed VPT ---------------------------- */

    // std::vector<double> X;

    // MPI_Init(NULL, NULL);

    // int world_size;
    // MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // int process_rank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
    
    // if (process_rank == 0) {
    //     X.resize(n * d);
    //     util::readToRowMajorVector(X, n, d, fileName);
    //     // prt::rowMajor(X.data(), n, d);
    // }

    // // std::cout << "process " << process_rank << "\t X.size() = " << X.size() << std::endl;

    // // Number of elements to distribute to each process
    // std::vector<int> chunk_size(world_size);
    // // The displacements where each chunk begins
    // std::vector<int> displs(world_size);

    // util::computeChunksDisplacements(chunk_size.data(), displs.data(), world_size, n, d);

    // // First chunk is the first to be increased in size in case world_size is not perfectly divisible by n
    // int MAX_CHUNK_S = chunk_size[0];
    // std::vector<double> _X(MAX_CHUNK_S);
    // std::vector<double> _Y(MAX_CHUNK_S);
    // std::vector<double> _Z(MAX_CHUNK_S);

    // /* ------------------------------ Distribute X ------------------------------ */

    // MPI_Scatterv(X.data(),
    //              chunk_size.data(),
    //              displs.data(),
    //              MPI_DOUBLE,
    //              _X.data(),
    //              chunk_size[process_rank],
    //              MPI_DOUBLE,
    //              0,
    //              MPI_COMM_WORLD);

    // // std::cout << "process " << process_rank << " displs = " << displs[process_rank] << std::endl;
    // // prt::rowMajor(_X.data(), chunk_size[process_rank] / d, d);

    // /* ----------------------------- Build local VPT ---------------------------- */

    // int* _indices = new int[chunk_size[process_rank] / d];
    // for (int i = 0; i < chunk_size[process_rank] / d; i++) {
    //     _indices[i] = i + displs[process_rank] / d;
    // }

    // std::vector<Point> _corpus(chunk_size[process_rank] / d);
    // conv::recVector(_corpus, _indices, _X.data(), d);
    // // prt::points(_corpus);
    
    // VPT _vpt(_corpus, b, k);
    // Node* _root = _vpt.buildTree(0, _corpus.size());
    // std::cout << "\nProcess " << process_rank << " local vantage point tree:\n\n";
    // prt::tree(_root, _corpus);
    
    // /* -------------------------------------------------------------------------- */




    // MPI_Finalize();

    /* -------------------------------------------------------------------------- */

    std::cout << std::endl;

    return 0;
}