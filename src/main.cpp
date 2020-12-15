#include <iostream>
#include <stdlib.h>

#include "knn.h"

int main() {

    std::cout << std::endl;
    
    const int n = 5;
    const int m = 3;
    const int d = 1;

    // double X[n * d] = {3, 3, 5, 1};
    // double Y[m * d] = {0, 1, 2, 5, 4, 1};

    double X[n * d] = {0, 10, -10, 20, 30};
    double Y[m * d] = {1, 3, 9};

    double D[m * n];

    const int k = 2;

    struct knnresult res = kNN(X, Y, D, n, m, d, k);

    std::cout << "kNN distances: " << std::endl;

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < k; j++)
            std::cout << res.ndist[i][j] << " ";
        std::cout << std::endl;
    }

    std::cout << std::endl << "kNN indices: " << std::endl;

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < k; j++)
            std::cout << res.nidx[i][j] << " ";
        std::cout << std::endl;
    }

    std::cout << std::endl;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++)
            std::cout << D[i * m + j] << " ";
        std::cout << std::endl;
    }

    std::cout << std::endl;
    
    return 0;
}
