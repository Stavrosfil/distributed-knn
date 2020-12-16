#ifndef __DISTRIBUTED_H__
#define __DISTRIBUTED_H__

#include <mpi.h>
#include <stdio.h>

#include "knn.hpp"

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

    std::cout << "test" << std::endl;

    return knnresult();
}

void test() {

  int matrix[] = {0, 1, 5, 2, 15, 7, 9150, 14};
  int N = 8;

  MPI_Init(NULL, NULL);

  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  int process_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);

  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);

  int *sub_matrix = new int[N/world_size];

  MPI_Scatter(matrix, N/world_size, MPI_INT, sub_matrix, N/world_size, MPI_INT, 0, MPI_COMM_WORLD);
  
  //sub_matrix[0] = 999;
  printf("%d %d ", sub_matrix[0], sub_matrix[1]);

  int buffer_len = 512;
  char buffer[buffer_len];
  sprintf(buffer, "Hello, Git! Rank: %d Total: %d Machine: %s", process_rank, world_size, processor_name);

  if (process_rank == 0) {

    //printf("%s\n", buffer);
    for (int i = 1; i < world_size; i++) {

      MPI_Recv(buffer, buffer_len, MPI_CHAR, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      //printf("%s\n", buffer);
    }
  }
  else {

    MPI_Send(buffer, buffer_len, MPI_CHAR, 0, process_rank, MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  //MPI_Gather(sub_matrix, N/world_size, MPI_INT, matrix, N/world_size, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Finalize();

  if (process_rank == 0) {

    // for (int i = 0; i < 8; i++) {

    // printf("%d ", matrix[i]);
    // }
  }

  //printf("\n");
  
  for (int i = 0; i < 8; i++) {

    //printf("%d ", matrix[i]);
  }

  //delete[] matrix;
  delete[] sub_matrix;

}

} // namespace mpi

#endif // __DISTRIBUTED_H__