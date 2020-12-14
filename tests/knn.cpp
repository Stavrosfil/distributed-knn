#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"
#include "knn.h"

// Tests that don't naturally fit in the headers/.cpp files directly
// can be placed in a tests/*.cpp file. Integration tests are a good example.

TEST_CASE("kNN with one dimensional points") {

    const int n = 7;
    const int m = 3;
    const int d = 1;
    const int k = 3;
    double D[n * m];

    SUBCASE("Distinct points bewteen X and Y") {

        double X[n * d] = {1.0, 2.0, 4.0, 15.0, 16.0, 18.0, 40.0};
        double Y[m * d] = {0.0, 17.0, 30.0};

        double ndist_ans[][k] = {1, 2, 4, 1, 1, 2, 10, 12, 14};
        int nidx_ans[][k]     = {0, 1, 2, 5, 4, 3, 6, 5, 4};

        struct knnresult res = kNN(X, Y, D, n, m, d, k);

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < k; j++) {
                CHECK(res.nidx[i][j] == nidx_ans[i][j]);
                CHECK(res.ndist[i][j] == ndist_ans[i][j]);
            }
        }
    }

    SUBCASE("One same point in X and Y") {
        double X[n * d] = {0.0, 2.0, 4.0, 15.0, 16.0, 18.0, 40.0};
        double Y[m * d] = {0.0, 17.0, 30.0};

        double ndist_ans[][k] = {2, 4, 15, 1, 1, 2, 10, 12, 14};
        int nidx_ans[][k]     = {1, 2, 3, 5, 4, 3, 6, 5, 4};

        struct knnresult res = kNN(X, Y, D, n, m, d, k);

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < k; j++) {
                CHECK(res.nidx[i][j] == nidx_ans[i][j]);
                CHECK(res.ndist[i][j] == ndist_ans[i][j]);
            }
        }
    }
}

TEST_CASE("kNN with two dimensional points") {

    const int n = 2;
    const int m = 3;
    const int d = 2;
    const int k = 2;
    double D[n * m];

    double X[n * d] = {3.0, 3.0, 5.0, 1.0};
    double Y[m * d] = {0.0, 1.0, 2.0, 5.0, 4.0, 1.0};

    double ndist_ans[][k] = {3.60555, 5, 2.23607, 5, 1, 2.23607};
    int nidx_ans[][k]     = {0, 1, 0, 1, 1, 0};

    struct knnresult res = kNN(X, Y, D, n, m, d, k);

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < k; j++) {
            CHECK(doctest::Approx(res.nidx[i][j]) == nidx_ans[i][j]);
            CHECK(doctest::Approx(res.ndist[i][j]) == ndist_ans[i][j]);
        }
    }
}
