#ifndef __DISTRIBUTED_H__
#define __DISTRIBUTED_H__

#include <mpi.h>
#include <stdio.h>

#include "knn.hpp"
#include "utils.hpp"

namespace mpi {

//! Compute distributed all-kNN of points in X
/*!

  \param  X      Data points                     [n-by-d]
  \param  n      Number of data points           [scalar]
  \param  d      Number of dimensions            [scalar]
  \param  k      Number of neighbors             [scalar]

  \return  The kNN result
*/
knnresult distrAllkNN(double* X, int n, int d, int k) {

    int N = n * d;

    MPI_Init(NULL, NULL);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int process_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);

    int chunk_size  = N / world_size;
    int xlen        = N / world_size;
    double* X_local = new double[xlen];

    MPI_Scatter(X, xlen, MPI_DOUBLE, X_local, xlen, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /* ------------------------------ Calculations ------------------------------ */

    // index = i + process_rank * chunk_size

    if (process_rank == world_size - 1)
        xlen = N - chunk_size * process_rank;

    for (int i = 0; i < xlen; i++) {
        std::cout << "Rank: " << process_rank << ", cs: " << xlen << " -> " << X_local[i] << std::endl;
    }

    int* buffer = new int[world_size];
    // struct knnresult* ans = new knnresult[world_size];
    struct knnresult ans = knnresult();

    for (int i = 0; i < world_size; i++) {
        // ans[i] = kNN(X_local, X_local, xlen, xlen, d, k);
        ans = kNN(X_local, X_local, i * chunk_size, xlen, xlen, d, k);

        MPI_Send(X_local, 1, MPI_DOUBLE, (process_rank + 1) % world_size, 0, MPI_COMM_WORLD);

        if (process_rank == 0) {
            MPI_Recv(X_local, 1, MPI_DOUBLE, world_size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            MPI_Recv(X_local, 1, MPI_DOUBLE, process_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    prt::twoDim(ans.nidx, xlen, k);
    prt::twoDim(ans.ndist, xlen, k);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    delete[] X_local;
    return knnresult();
}

} // namespace mpi

#endif // __DISTRIBUTED_H__