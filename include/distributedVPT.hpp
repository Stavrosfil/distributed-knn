#ifndef __DISTRIBUTEDVPT_H__
#define __DISTRIBUTEDVPT_H__

#include <mpi.h>
#include <stdio.h>

#include "knn.hpp"
#include "utils.hpp"

namespace mpi {

knnresult distrVPTkNN(std::vector<double> X, int n, int d, int k, int b, std::string fileName)
{

    util::Timer timer(false);

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
    std::vector<int> _Xindices(MAX_CHUNK_S / d);
    std::vector<double> _Xcoords(MAX_CHUNK_S);
    std::vector<int> _Yindices(chunk_size[process_rank] / d);
    std::vector<double> _Ycoords(chunk_size[process_rank]);
    std::vector<int> _Zindices(MAX_CHUNK_S / d);
    std::vector<double> _Zcoords(MAX_CHUNK_S);

    /* ------------------------------ Distribute X ------------------------------ */

    // cl::start();

    MPI_Scatterv(X.data(),
                 chunk_size.data(),
                 displs.data(),
                 MPI_DOUBLE,
                 _Xcoords.data(),
                 chunk_size[process_rank],
                 MPI_DOUBLE,
                 0,
                 MPI_COMM_WORLD);

    // cl::stop(true, "Data distribution");
    // std::cout << "process " << process_rank << " displs = " << displs[process_rank] << std::endl;
    // prt::rowMajor(_Xcoords.data(), MAX_CHUNK_S / d, d);

    /* ----------------------------- Build local VPT ---------------------------- */

    for (int i = 0; i < _Xindices.size(); i++) {
        _Xindices[i] = i + displs[process_rank] / d;
    }

    std::vector<Point> _corpus(MAX_CHUNK_S / d);
    conv::recVector(_corpus, _Xindices, _Xcoords, d, chunk_size[process_rank] / d);

    // _query is prefered to have the initial sort (before build tree) for simplicity in kNN
    std::vector<Point> _query(_corpus);

    // prt::points(_query);

    VPT _vpt(_corpus, b, k, chunk_size[process_rank] / d);
    _vpt.build(chunk_size[process_rank] / d);

    // we need the sorted form of _X for tree reconstruction
    conv::serVector(_corpus, _Xindices, _Xcoords, chunk_size[process_rank] / d);

    // prt::points(_corpus);
    // prt::points(_query);

    // std::cout << "\nProcess " << process_rank << " local vantage point tree:\n\n";
    // prt::tree(_root, _corpus);

    /* ------------------------------ Calculations ------------------------------ */

    knnresult _ans = knnresult();
    _ans.m         = chunk_size[process_rank] / d;
    _ans.k         = k;
    std::vector<int> nidx(_ans.m * _ans.k, -1);
    std::vector<double> ndist(_ans.m * _ans.k, D_MAX);
    _ans.nidx  = nidx.data();
    _ans.ndist = ndist.data();

    int prev_rank              = (world_size + process_rank - 1) % world_size;
    int next_rank              = (world_size + process_rank + 1) % world_size;
    unsigned rolling_prev_rank = world_size + process_rank + 1;

    timer.start("MPI Ring");

    for (int i = 0; i < world_size; i++) {

        rolling_prev_rank = (rolling_prev_rank + world_size - 1) % world_size;

        MPI_Request* reqs = new MPI_Request[4];

        MPI_Irecv(/* recv buffer: */ _Zindices.data(),
                  /* count: */ _Zindices.size(),
                  /* type: */ MPI_INT,
                  /* from: */ prev_rank,
                  /* tag: */ 0,
                  /* communicator: */ MPI_COMM_WORLD,
                  /* request: */ &reqs[0]);
        MPI_Irecv(/* recv buffer: */ _Zcoords.data(),
                  /* count: */ _Zcoords.size(),
                  /* type: */ MPI_DOUBLE,
                  /* from: */ prev_rank,
                  /* tag: */ 0,
                  /* communicator: */ MPI_COMM_WORLD,
                  /* request: */ &reqs[1]);

        MPI_Isend(/* send buffer: */ _Xindices.data(),
                  /* count: */ _Xindices.size(),
                  /* type: */ MPI_INT,
                  /* to: */ next_rank,
                  /* tag: */ 0,
                  /* communicator: */ MPI_COMM_WORLD,
                  /* request: */ &reqs[2]);
        MPI_Isend(/* send buffer: */ _Xcoords.data(),
                  /* count: */ _Xcoords.size(),
                  /* type: */ MPI_DOUBLE,
                  /* to: */ next_rank,
                  /* tag: */ 0,
                  /* communicator: */ MPI_COMM_WORLD,
                  /* request: */ &reqs[3]);

        // reconstruct vector
        // update _corpus with the received serialized _X
        conv::recVector(_corpus, _Xindices, _Xcoords, d, chunk_size[rolling_prev_rank] / d);

        _vpt.rebuild(chunk_size[rolling_prev_rank] / d);

        // std::cout << "Process rank ->\t" << process_rank << std::endl;
        // if (process_rank == 1)
        // prt::points(_corpus);

        // if (process_rank) {
        //     std::cout << "\nprocess " << process_rank << " reconstructed vpt:\n";
        //     _vpt.printTree();
        // }

        // kNN
        for (int j = 0; j < chunk_size[process_rank] / d; j++) {
            _vpt.computeKNN(_query[j], _ans, j);
        }

        // std::cout << "process " << process_rank << " local kNN:\n";
        // prt::kNN(_ans);

        MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);

        _Xindices = _Zindices;
        _Xcoords  = _Zcoords;
    }

    timer.stop();

    /* ----------------------------- Gather results ----------------------------- */

    knnresult ans = knnresult();

    if (process_rank == 0) {
        ans.m     = n;
        ans.k     = k;
        ans.nidx  = new int[n * k];
        ans.ndist = new double[n * k];
    }

    int* recv_chunks = new int[world_size];
    int* recv_displs = new int[world_size];

    util::computeChunksDisplacements(recv_chunks, recv_displs, world_size, n, k);

    MPI_Gatherv(_ans.nidx, _ans.m * _ans.k, MPI_INT, ans.nidx, recv_chunks, recv_displs, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gatherv(
        _ans.ndist, _ans.m * _ans.k, MPI_DOUBLE, ans.ndist, recv_chunks, recv_displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // if (process_rank == 0)
    // prt::kNN(ans);

    /* -------------------------------------------------------------------------- */

    MPI_Finalize();

    return ans;

    MPI_Finalize();
    return knnresult();
}

} // namespace mpi

#endif // __DISTRIBUTEDVPT_H__