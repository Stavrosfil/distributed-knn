#ifndef __DISTRIBUTEDVPT_H__
#define __DISTRIBUTEDVPT_H__

#include <mpi.h>
#include <stdio.h>

#include "knn.hpp"
#include "utils.hpp"

namespace mpi {

void distrVPTkNN(std::vector<double> X, int n, int d, int k, int b, std::string fileName)         // TODO knnresult return
{
    /* --------------------------- Init Communication --------------------------- */

    MPI_Init(NULL, NULL);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int process_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
    
    if (process_rank == 0) {
        X.resize(n * d);
        util::readToRowMajorVector(X, n, d, fileName);
        // prt::rowMajor(X.data(), n, d);
    }

    // std::cout << "process " << process_rank << "\t X.size() = " << X.size() << std::endl;

    // Number of elements to distribute to each process
    std::vector<int> chunk_size(world_size);
    // The displacements where each chunk begins
    std::vector<int> displs(world_size);

    util::computeChunksDisplacements(chunk_size.data(), displs.data(), world_size, n, d);

    // First chunk is the first to be increased in size in case world_size is not perfectly divisible by n
    int MAX_CHUNK_S = chunk_size[0];
    std::vector<double> _X(MAX_CHUNK_S);
    std::vector<double> _Y(MAX_CHUNK_S);
    std::vector<double> _Z(MAX_CHUNK_S);

    /* ------------------------------ Distribute X ------------------------------ */

    // cl::start();

    MPI_Scatterv(X.data(),
                 chunk_size.data(),
                 displs.data(),
                 MPI_DOUBLE,
                 _X.data(),
                 chunk_size[process_rank],
                 MPI_DOUBLE,
                 0,
                 MPI_COMM_WORLD);

    // cl::stop(true, "Data distribution");
    // std::cout << "process " << process_rank << " displs = " << displs[process_rank] << std::endl;
    // prt::rowMajor(_X.data(), chunk_size[process_rank] / d, d);

    /* ----------------------------- Build local VPT ---------------------------- */

    // cl::start();

    int* _indices = new int[chunk_size[process_rank] / d];
    for (int i = 0; i < chunk_size[process_rank] / d; i++) {
        _indices[i] = i + displs[process_rank] / d;
    }

    std::vector<Point> _corpus(chunk_size[process_rank] / d);
    conv::recVector(_corpus, _indices, _X.data(), d);
    // prt::points(_corpus);
    
    VPT _vpt(_corpus, b, k);
    Node* _root = _vpt.buildTree(0, _corpus.size());

    // cl::stop(true, "Build local tree");
    // std::cout << "\nProcess " << process_rank << " local vantage point tree:\n\n";
    // prt::tree(&_root, _corpus);
    
    /* ------------------------------ Calculations ------------------------------ */

    _Y = _X;

    std::vector<Point> _query(_corpus);

    // prt::points(_query);

    knnresult _ans = knnresult();
    _ans.m         = chunk_size[process_rank] / d;
    _ans.k         = k;
    std::vector<int> nidx(_ans.m * _ans.k, -1);
    std::vector<double> ndist(_ans.m * _ans.k, D_MAX);
    _ans.nidx = nidx.data();
    _ans.ndist = ndist.data();

    int prev_rank              = (world_size + process_rank - 1) % world_size;
    int next_rank              = (world_size + process_rank + 1) % world_size;
    unsigned rolling_prev_rank = world_size + process_rank + 1;

    for (auto p : _query) {
        _vpt.kNN(p, _ans, p.index - displs[process_rank] / d, *_root);
    }

    // prt::kNN(_ans);

    /* -------------------------------------------------------------------------- */

    MPI_Finalize();

}

} // namespace mpi

#endif // __DISTRIBUTEDVPT_H__