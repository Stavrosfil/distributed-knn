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

    int rem         = n % world_size;      // elements remaining after division among processes
    int sum         = 0;                   // Sum of counts. Used to calculate displacements
    int* chunk_size = new int[world_size]; // array describing how many elements to send to each process
    int* displs     = new int[world_size]; // array describing the displacements where each segment begins

    // calculate chunk sizes and displacements
    for (int i = 0; i < world_size; i++) {
        chunk_size[i] = n / world_size;
        if (rem > 0) {
            chunk_size[i]++;
            rem--;
        }
        displs[i] = sum;
        sum += chunk_size[i];
    }

    // First chunk is the first to be increased in size in case world_size is not perfectly divisible by n
    int MAX_CHUNK_S = chunk_size[0];

    double* _X = new double[chunk_size[process_rank]]; // chunk_size[process_rank] = n/world_size OR n/world_size+1
    double* _Y = new double[MAX_CHUNK_S];              // set maximum length
    double* _Z = new double[MAX_CHUNK_S];

    /* ------------------------------ Distribute X ------------------------------ */

    MPI_Scatterv(X, chunk_size, displs, MPI_DOUBLE, _X, chunk_size[process_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /* ------------------------------ Calculations ------------------------------ */

    _Y                    = _X;
    struct knnresult _res = knnresult();
    _res.m                = chunk_size[process_rank];
    _res.k                = k;
    _res.nidx             = new int[_res.m * _res.k]();
    _res.ndist            = new double[_res.m * _res.k]();

    for (int i = 0; i < _res.m; i++) {
        for (int j = 0; j < k; j++) {
            _res.ndist[i * k + j] = D_MAX;
            _res.nidx[i * k + j]  = -1;
        }
    }

    // For each process we iterate in a circular fashion
    for (int i = process_rank; i < process_rank + world_size; i++) {

        int prev_rank = (world_size + i - 1) % world_size;
        int next_rank = (i + 1) % world_size;

        // std::cout << "Process " << process_rank << " -> " << i << " " << next_rank << " " << prev_rank << std::endl;

        kNN(_res, _Y, _X, displs[i], chunk_size[i % world_size], chunk_size[process_rank], d, k);

        MPI_Send(_Y, MAX_CHUNK_S, MPI_DOUBLE, next_rank, 1, MPI_COMM_WORLD);
        MPI_Recv(_Z, MAX_CHUNK_S, MPI_DOUBLE, prev_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        _Y = _Z;
    }

    std::cout << "Process " << process_rank << std::endl;
    prt::twoDim(_res.ndist, _res.m, _res.k);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    // delete[] _X;
    // delete[] _Y;
    // delete[] _Z;
    return _res;
}

} // namespace mpi

#endif // __DISTRIBUTED_H__