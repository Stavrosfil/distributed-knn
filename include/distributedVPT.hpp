#ifndef __DISTRIBUTEDVPT_H__
#define __DISTRIBUTEDVPT_H__

#include <mpi.h>
#include <stdio.h>

#include "knn.hpp"
#include "utils.hpp"

namespace mpi {

knnresult distrVPTkNN(std::vector<double> X, int n, int d, int k, int b, std::string fileName)
{
    util::Timer timer(true);

    /* --------------------------- Init Communication --------------------------- */

    MPI_Init(NULL, NULL);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int process_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);

    if (process_rank == 0) {
        std::cout << std::endl;
        X.resize(n * d);
        util::read(X, n, d, fileName);
        // prt::rowMajor(X.data(), n, d);
        timer.start("v2");
    }

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

    MPI_Scatterv(X.data(),
                 chunk_size.data(),
                 displs.data(),
                 MPI_DOUBLE,
                 _Xcoords.data(),
                 chunk_size[process_rank],
                 MPI_DOUBLE,
                 0,
                 MPI_COMM_WORLD);

    /* ----------------------------- Build local VPT ---------------------------- */

    for (int i = 0; i < _Xindices.size(); i++) {
        _Xindices[i] = i + displs[process_rank] / d;
    }

    std::vector<Point> _corpus(MAX_CHUNK_S / d);
    conv::recVector(_corpus, _Xindices, _Xcoords, d, chunk_size[process_rank] / d);

    // _query is prefered to have the initial sort (before build tree) for simplicity in kNN
    std::vector<Point> _query(_corpus);

    VPT _vpt(_corpus, b, k, chunk_size[process_rank] / d);
    _vpt.build(chunk_size[process_rank] / d);

    // we need the sorted form of _X for tree reconstruction
    conv::serVector(_corpus, _Xindices, _Xcoords, chunk_size[process_rank] / d);

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

    double nodes_visits = 0;
    // b = _Yindices.size() / 8500;

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
                  /* tag: */ 1,
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
                  /* tag: */ 1,
                  /* communicator: */ MPI_COMM_WORLD,
                  /* request: */ &reqs[3]);

        // reconstruct vector
        // update _corpus with the received serialized _X
        conv::recVector(_corpus, _Xindices, _Xcoords, d, chunk_size[rolling_prev_rank] / d);

        _vpt.rebuild(chunk_size[rolling_prev_rank] / d);

        // kNN
        for (int j = 0; j < chunk_size[process_rank] / d; j++) {
            _vpt.computeKNN(_query[j], _ans, j);
            nodes_visits += _vpt.getNodesVisits();
        }

        MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);

        _Xindices = _Zindices;
        _Xcoords  = _Zcoords;
    }

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

    /* -------------------------------------------------------------------------- */

    if (process_rank == 0) {
        timer.stop();
        // prt::kNN(ans);
    }

    nodes_visits /= _Yindices.size();
    
    // send all local nodes_visits to p0 and compute average
    if (process_rank != 0) {
        MPI_Send( /* data: */ &nodes_visits,
                  /* count: */ 1,
                  /* type: */ MPI_DOUBLE,
                  /* to: */ 0,
                  /* tag: */ process_rank,
                  /* communicator: */ MPI_COMM_WORLD);
    }
    else {
        std::vector<double> all_nodes_visits(world_size);
        all_nodes_visits[0] = nodes_visits;
        
        for (int i = 1; i < world_size; i ++) {
            MPI_Recv( /* data: */ &all_nodes_visits[i],
                      /* count: */ 1,
                      /* type: */ MPI_DOUBLE,
                      /* from: */ i,
                      /* tag: */ i,
                      /* communicator: */ MPI_COMM_WORLD,
                      /* status: */ MPI_STATUS_IGNORE);
        }

        double av_nodes_visits = 0;

        for (auto& n : all_nodes_visits)
            av_nodes_visits += n;
        
        av_nodes_visits /= world_size;

        std::cout << "average nodes visits per query point = " << av_nodes_visits << std::endl;
    }

    if (process_rank == 0)
        std::cout << std::endl;

    MPI_Finalize();

    return ans;
}

} // namespace mpi

#endif // __DISTRIBUTEDVPT_H__