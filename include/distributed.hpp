#ifndef __DISTRIBUTED_H__
#define __DISTRIBUTED_H__

#include <mpi.h>
#include <stdio.h>

#include "knn.hpp"
#include "utils.hpp"

namespace mpi {

knnresult distrAllkNN(std::vector<double> X, int n, int d, int k, std::string fileName)
{

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

    // Number of elements to distribute to each process
    int* chunk_size = new int[world_size];
    // The displacements where each chunk begins
    int* displs = new int[world_size];

    util::computeChunksDisplacements(chunk_size, displs, world_size, n, d);

    // First chunk is the first to be increased in size in case world_size is not perfectly divisible by n
    int MAX_CHUNK_S = chunk_size[0];
    double* _X      = new double[chunk_size[process_rank]];
    double* _Y      = new double[MAX_CHUNK_S];
    double* _Z      = new double[MAX_CHUNK_S];

    /* ------------------------------ Distribute X ------------------------------ */

    MPI_Scatterv(X.data(), chunk_size, displs, MPI_DOUBLE, _X, chunk_size[process_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /* ------------------------------ Calculations ------------------------------ */
    
    memcpy(_Y, _X, MAX_CHUNK_S * sizeof(double));

    knnresult _ans = knnresult();
    _ans.m         = chunk_size[process_rank] / d;
    _ans.k         = k;
    _ans.nidx      = new int[_ans.m * _ans.k]();
    _ans.ndist     = new double[_ans.m * _ans.k]();

    for (int i = 0; i < _ans.m; i++) {
        for (int j = 0; j < k; j++) {
            _ans.ndist[i * k + j] = D_MAX;
            _ans.nidx[i * k + j]  = -1;
        }
    }

    int prev_rank              = (world_size + process_rank - 1) % world_size;
    int next_rank              = (world_size + process_rank + 1) % world_size;
    unsigned rolling_prev_rank = world_size + process_rank + 1;

    for (int i = 0; i < world_size; i++) {

        rolling_prev_rank = --rolling_prev_rank % world_size;

        MPI_Request* reqs = new MPI_Request[2];

        MPI_Irecv(/* recv buffer: */ _Z,
                  /* count: */ MAX_CHUNK_S,
                  /* type: */ MPI_DOUBLE,
                  /* from: */ prev_rank,
                  /* tag: */ 0,
                  /* communicator: */ MPI_COMM_WORLD,
                  /* request: */ &reqs[0]);
        MPI_Isend(/* send buffer: */ _Y,
                  /* count: */ MAX_CHUNK_S,
                  /* type: */ MPI_DOUBLE,
                  /* to: */ next_rank,
                  /* tag: */ 0,
                  /* communicator: */ MPI_COMM_WORLD,
                  /* request: */ &reqs[1]);

        kNN(_ans,
            _Y,
            _X,
            displs[rolling_prev_rank] / d,
            chunk_size[rolling_prev_rank] / d,
            chunk_size[process_rank] / d,
            d,
            k);

        MPI_Waitall(2, reqs, MPI_STATUSES_IGNORE);

        memcpy(_Y, _Z, MAX_CHUNK_S * sizeof(double));
    }

    knnresult ans = knnresult();
    ans.m         = n;
    ans.k         = k;
    ans.nidx      = new int[n * k];
    ans.ndist     = new double[n * k];

    int* recv_chunks = new int[world_size];
    int* recv_displs = new int[world_size];

    util::computeChunksDisplacements(recv_chunks, recv_displs, world_size, n, k);

    MPI_Gatherv(_ans.nidx, _ans.m * _ans.k, MPI_INT, ans.nidx, recv_chunks, recv_displs, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gatherv(
        _ans.ndist, _ans.m * _ans.k, MPI_DOUBLE, ans.ndist, recv_chunks, recv_displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Finalize();

    if (process_rank == 0)
        prt::kNN(ans);
        
    // if (process_rank == 0)
    //     prt::kNN(_ans);

    return ans;
}

} // namespace mpi

#endif // __DISTRIBUTED_H__