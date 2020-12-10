#include <iostream>  
#include <stdio.h>     
#include <cmath> 
#include <limits>

using namespace std; 


// Definition of the kNN result struct
typedef struct knnresult{

    int    * nidx;    //!< Indices (0-based) of nearest neighbors [m-by-k]
    double * ndist;   //!< Distance of nearest neighbors          [m-by-k]
    int      m;       //!< Number of query points                 [scalar]
    int      k;       //!< Number of nearest neighbors            [scalar]
    
} knnresult;


// add function
// Add an element to knn arrays and sort array.
void add(double dist, double* ndist, int* nidx, int m, int Xi, int k) {

        for (int i = k - 1; i > 0; i--) {
            if (ndist[i - 1] < dist) {
                for (int j = k - 1; j > i; j--) {
                    ndist[j] = ndist[j - 1];
                    nidx[j] = nidx[j - 1];
                }
                ndist[i] = dist;
                nidx[i] = Xi;
                break;
            }
        }
        if (ndist[0] >= dist) {
            for (int j = k - 1; j > 0; j--) {
                    ndist[j] = ndist[j - 1];
                    nidx[j] = nidx[j - 1];
            }
            ndist[0] = dist;
            nidx[0] = Xi;
        }

    }



int main() {

    // int n = 2;
    // int m = 3;
    // int d = 2;
    // int k = 2;
    // int X[n][d] = {{3, 3}, {5, 1}};
    // int Y[m][d] = {{0, 1}, {2, 5}, {4, 1}};

    int n = 7;
    int m = 3;
    int d = 1;
    int k = 3;
    int X[n][d] = {{1}, {2}, {4}, {15}, {16}, {18}, {40}};
    int Y[n][d] = {{0}, {17}, {30}};

    int nidx[m][k];
    double ndist[m][k];

    // init ndist
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < k; j++) {
            ndist[i][j] = std::numeric_limits<double>::max();
        }
    }

    // find knn
    for (int i = 0; i < m; i++) {                                   // query data points Y
        for (int j = 0; j < n; j++) {                               // corpus data points X
            // int temp1 = Y[i][0] - X[j][0];
            // int temp2 = Y[i][1] - X[j][1];
            // double dist = sqrt(temp1*temp1 + temp2*temp2);          // calculate Yi - Xj distance for d = 2
            double dist = abs(Y[i][0] - X[j][0]);                      // calculate Yi - Xj distance for d = 1
            //cout << dist[i][j] << " ";
            if (dist < ndist[i][k - 1]) {
                add(dist, ndist[i], nidx[i], m, j, k);                 // add Xj and short knn arrays
            }
        } 
    }

    // print results
    cout << "knn distances: " << endl;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < k; j++) {
            cout << ndist[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl << "knn indeces: " << endl;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < k; j++) {
            cout << nidx[i][j] << " ";
        }
        cout << endl;
    }

    //knnresult kNN(double * X, double * Y, int n, int m, int d, int k);

    return 0;

}