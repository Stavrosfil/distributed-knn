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

    int MAX_CHUNK_S = n / world_size + 1;

    double* _X = new double[chunk_size[process_rank]]; // chunk_size[process_rank] = n/world_size OR n/world_size+1
    double* _Y = new double[MAX_CHUNK_S];              // set maximum length
    double* _Z = new double[MAX_CHUNK_S];

    MPI_Scatterv(X, chunk_size, displs, MPI_DOUBLE, _X, chunk_size[process_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    printf("Process %d _X: ", process_rank);

    for (int i = 0; i < chunk_size[process_rank]; i++) {
        printf("%.3lf\t", _X[i]);
    }
    std::cout << std::endl;

    /* ------------------------------ Calculations ------------------------------ */

    // for (int i = 0; i < xlen; i++) {
    //     std::cout << "Rank: " << process_rank << ", cs: " << xlen << " -> " << _X[i] << std::endl;
    // }

    _Y                    = _X;
    struct knnresult _res = knnresult();
    _res.m                = chunk_size[process_rank];
    _res.k                = k;
    _res.nidx             = new int[_res.m * _res.k];
    _res.ndist            = new double[_res.m * _res.k];

    // TODO : Create mpi data type "smessage" and "rmessage"
    // TODO : Create a struct "smessage" that will consist of: "struct knnresult distr_ans" and  "double * _Y"
    // TODO : Create a struct "rmessage" that will consist of: "struct knnresult distr_ans" and  "double * _Z"
    // NOTE : "struct knnresult distr_ans" includes the chunk_size
    // NOTE : The communication will be achived by sending and receiving just one struct

    for (int i = 0; i < _res.m; i++) {
        for (int j = 0; j < k; j++) {
            _res.ndist[i * k + j] = D_MAX;
        }
    }

    int offset = 0;

    // For each process we iterate in a circular fashion
    for (int i = process_rank; i < process_rank + world_size; i++) {

        std::cout << "********** Iteration " << i << " ************" << std::endl;

        kNN(_res, _X, _Y, offset, chunk_size[process_rank], MAX_CHUNK_S, d, k);

        MPI_Send(_Y, MAX_CHUNK_S, MPI_DOUBLE, (i + 1) % world_size, 1, MPI_COMM_WORLD);
        // MPI_Send(distr_ans, 1, mpi_knnresult_type, (process_rank + 1) % world_size, 02, MPI_COMM_WORLD);

        MPI_Recv(_Z, MAX_CHUNK_S, MPI_DOUBLE, (world_size + i - 1) % world_size, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // MPI_Recv(distr_ans, 1, mpi_knnresult_type, (process_rank + 1) % world_size, 02, MPI_COMM_WORLD);

        _Y = _Z;

        std::cout << "Process " << process_rank << " _Y updated: " << std::endl;

        for (int i = 0; i < MAX_CHUNK_S; i++) {

            // TODO : After the creation of the struct, "MAX_CHUNK_S" will
            // be replaced by "rmessage.distr_ans.m"

            std::cout << _Y[i] << " ";
        }
        std::cout << std::endl;
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ STAVROS OLD VERSION
    // FOLLOWING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // // struct knnresult* ans = new knnresult[world_size];
    // struct knnresult ans = knnresult();

    // for (int i = 0; i < world_size; i++) {
    //     // ans[i] = kNN(_X, _X, xlen, xlen, d, k);
    //     ans = kNN(_X, _X, i * chunk_size, xlen, xlen, d, k);

    //     MPI_Send(_X, 1, MPI_DOUBLE, (process_rank + 1) % world_size, 0, MPI_COMM_WORLD);

    //     if (process_rank == 0) {
    //         MPI_Recv(_X, 1, MPI_DOUBLE, world_size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //     } else {
    //         MPI_Recv(_X, 1, MPI_DOUBLE, process_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //     }
    // }

    // prt::twoDim(ans.nidx, xlen, k);
    // prt::twoDim(ans.ndist, xlen, k);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    // delete[] _X;
    // delete[] _Y;
    // delete[] _Z;
    return knnresult();
}

} // namespace mpi

#endif // __DISTRIBUTED_H__