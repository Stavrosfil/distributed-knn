#ifndef __DISTRIBUTED_H__
#define __DISTRIBUTED_H__

#include <mpi.h>
#include <stdio.h>

#include "knn.hpp"
#include "utils.hpp"

namespace mpi {

knnresult distrAllkNN(double* X, int n, int d, int k) {

    MPI_Init(NULL, NULL);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int process_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);

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

    MPI_Scatterv(X, chunk_size, displs, MPI_DOUBLE, _X, chunk_size[process_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /* ------------------------------ Calculations ------------------------------ */

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

    memcpy(_Y, _X, MAX_CHUNK_S * sizeof(double));
    kNN(_ans, _Y, _X, displs[process_rank] / d, chunk_size[process_rank] / d, chunk_size[process_rank] / d, d, k);

    int prev_rank              = (world_size + process_rank - 1) % world_size;
    int next_rank              = (world_size + process_rank + 1) % world_size;
    unsigned rolling_prev_rank = world_size + process_rank;

    for (int i = 0; i < world_size - 1; i++) {

        rolling_prev_rank = --rolling_prev_rank % world_size;

        if (process_rank == 0) {
            MPI_Recv(_Z, MAX_CHUNK_S, MPI_DOUBLE, prev_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(_Y, MAX_CHUNK_S, MPI_DOUBLE, next_rank, 1, MPI_COMM_WORLD);
        } else {
            MPI_Send(_Y, MAX_CHUNK_S, MPI_DOUBLE, next_rank, 1, MPI_COMM_WORLD);
            MPI_Recv(_Z, MAX_CHUNK_S, MPI_DOUBLE, prev_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        memcpy(_Y, _Z, MAX_CHUNK_S * sizeof(double));

        kNN(_ans, _Y, _X, displs[rolling_prev_rank] / d, chunk_size[rolling_prev_rank] / d,
            chunk_size[process_rank] / d, d, k);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    knnresult ans = knnresult();
    ans.m         = n;
    ans.k         = k;
    ans.nidx      = new int[n * k];
    ans.ndist     = new double[n * k];

    int* recv_chunks = new int[world_size];
    int* recv_displs = new int[world_size];

    util::computeChunksDisplacements(recv_chunks, recv_displs, world_size, n, k);

    MPI_Gatherv(_ans.nidx, _ans.m * _ans.k, MPI_INT, ans.nidx, recv_chunks, recv_displs, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gatherv(_ans.ndist, _ans.m * _ans.k, MPI_DOUBLE, ans.ndist, recv_chunks, recv_displs, MPI_DOUBLE, 0,
                MPI_COMM_WORLD);

    MPI_Finalize();

    if (process_rank == 0)
        prt::rowMajor(ans.ndist, ans.m, k);

    return ans;
}

} // namespace mpi

#endif // __DISTRIBUTED_H__