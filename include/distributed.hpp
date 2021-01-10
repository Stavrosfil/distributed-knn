#ifndef __DISTRIBUTED_H__
#define __DISTRIBUTED_H__

#include <mpi.h>
#include <stdio.h>

#include "knn.hpp"
#include "utils.hpp"
#include "reader.hpp"

namespace mpi {

knnresult distrAllkNN(std::vector<double> X, int n, int d, int k, std::string fileName)
{
    util::Timer timer(true);

    /* --------------------------- Init Communication --------------------------- */
    MPI_Init(NULL, NULL);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int process_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);

    rdCorel::colorMom(n, d, X, process_rank);

    if (process_rank == 0)
        timer.start("v1");
    
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

    MPI_Scatterv(X.data(),
                chunk_size.data(),
                displs.data(),
                MPI_DOUBLE,
                _X.data(),
                chunk_size[process_rank],
                MPI_DOUBLE,
                0,
                MPI_COMM_WORLD);

/* ------------------------------ Calculations ------------------------------ */

    _Y = _X;

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

    for (int i = 0; i < world_size; i++) {

        rolling_prev_rank = (rolling_prev_rank + world_size - 1) % world_size;

        MPI_Request* reqs = new MPI_Request[2];

        MPI_Irecv(/* recv buffer: */ _Z.data(),
                /* count: */ MAX_CHUNK_S,
                /* type: */ MPI_DOUBLE,
                /* from: */ prev_rank,
                /* tag: */ 0,
                /* communicator: */ MPI_COMM_WORLD,
                /* request: */ &reqs[0]);
        MPI_Isend(/* send buffer: */ _Y.data(),
                /* count: */ MAX_CHUNK_S,
                /* type: */ MPI_DOUBLE,
                /* to: */ next_rank,
                /* tag: */ 0,
                /* communicator: */ MPI_COMM_WORLD,
                /* request: */ &reqs[1]);

        kNN(_ans,
            _Y.data(),
            _X.data(),
            displs[rolling_prev_rank] / d,
            chunk_size[rolling_prev_rank] / d,
            chunk_size[process_rank] / d,
            d,
            k);

        MPI_Waitall(2, reqs, MPI_STATUSES_IGNORE);

        _Y = _Z;
    }

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

    if (process_rank == 0) {
        timer.stop();
        std::cout << std::endl;
        //prt::kNN(ans);
    }

    MPI_Finalize();

    return ans;
}

} // namespace mpi

#endif // __DISTRIBUTED_H__