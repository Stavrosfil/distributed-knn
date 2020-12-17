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

    int rem = n % world_size;               // elements remaining after division among processes
    int sum = 0;                            // Sum of counts. Used to calculate displacements
    int *chunk_size = new int[world_size];  // array describing how many elements to send to each process
    int *displs = new int[world_size];      // array describing the displacements where each segment begins

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

    double *X_local = new double[chunk_size[process_rank]];         // chunk_size[process_rank] = n/world_size OR n/world_size+1
    double *Y_local = new double[n / world_size + 1];               // set maximum length
    double *Z_local = new double[n / world_size + 1];
    
    MPI_Scatterv(X, chunk_size, displs, MPI_DOUBLE, X_local, chunk_size[process_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    // print chunks
    printf("Process %d X_local: ", process_rank);
    for (int i = 0; i < chunk_size[process_rank]; i++) {
        printf("%lf\t", X_local[i]);
    }
    printf("\n");

    /* ------------------------------ Calculations ------------------------------ */

    // index = i + process_rank * chunk_size

    // for (int i = 0; i < xlen; i++) {
    //     std::cout << "Rank: " << process_rank << ", cs: " << xlen << " -> " << X_local[i] << std::endl;
    // }

    Y_local = X_local;
    struct knnresult distr_ans = knnresult();
    distr_ans.m = chunk_size[process_rank];
    distr_ans.k = k;

// TODO : Create mpi data type "smessage" and "rmessage"
// TODO : Create a struct "smessage" that will consist of: "struct knnresult distr_ans" and  "double * Y_local" 
// TODO : Create a struct "rmessage" that will consist of: "struct knnresult distr_ans" and  "double * Z_local" 
// NOTE : "struct knnresult distr_ans" includes the chunk_size 
// NOTE : The communication will be achived by sending and receiving just one struct

    MPI_Send(Y_local, n / world_size + 1, MPI_DOUBLE, (process_rank + 1) % world_size, 01, MPI_COMM_WORLD);
    //MPI_Send(distr_ans, 1, mpi_knnresult_type, (process_rank + 1) % world_size, 02, MPI_COMM_WORLD);

    MPI_Recv(Z_local, n / world_size + 1, MPI_DOUBLE, (process_rank + 1) % world_size, 01, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //MPI_Recv(distr_ans, 1, mpi_knnresult_type, (process_rank + 1) % world_size, 02, MPI_COMM_WORLD);

    Y_local = Z_local;

    std::cout << "Process " << process_rank << " Y_local updated: " << std::endl;
    for (int i = 0; i < n / world_size + 1; i++) {          // TODO : After the creation of the struct, "n / world_size + 1" will be replaced by "rmessage.distr_ans.m"
        
        std::cout << Y_local[i] << " ";
    }
    std::cout << std::endl;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ STAVROS OLD VERSION FOLLOWING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // // struct knnresult* ans = new knnresult[world_size];
    // struct knnresult ans = knnresult();

    // for (int i = 0; i < world_size; i++) {
    //     // ans[i] = kNN(X_local, X_local, xlen, xlen, d, k);
    //     ans = kNN(X_local, X_local, i * chunk_size, xlen, xlen, d, k);

    //     MPI_Send(X_local, 1, MPI_DOUBLE, (process_rank + 1) % world_size, 0, MPI_COMM_WORLD);

    //     if (process_rank == 0) {
    //         MPI_Recv(X_local, 1, MPI_DOUBLE, world_size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //     } else {
    //         MPI_Recv(X_local, 1, MPI_DOUBLE, process_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //     }
    // }

    // prt::twoDim(ans.nidx, xlen, k);
    // prt::twoDim(ans.ndist, xlen, k);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    // delete[] X_local;
    // delete[] Y_local;
    // delete[] Z_local;
    return knnresult();
}

} // namespace mpi

#endif // __DISTRIBUTED_H__