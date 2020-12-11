#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"
#include "knn.h"

// Tests that don't naturally fit in the headers/.cpp files directly
// can be placed in a tests/*.cpp file. Integration tests are a good example.

TEST_CASE("kNN with one dimensional points") {

    int n = 7;
    int m = 3;
    int d = 1;
    int k = 3;

    SUBCASE("Distinct points bewteen X and Y") {

        double X[][1] = {{1.0}, {2.0}, {4.0}, {15.0}, {16.0}, {18.0}, {40.0}};
        double Y[][1] = {{0.0}, {17.0}, {30.0}};

        double ndist_ans[][3] = {{1, 2, 4}, {1, 1, 2}, {10, 12, 14}};
        int nidx_ans[][3]     = {{0, 1, 2}, {5, 4, 3}, {6, 5, 4}};

        struct knnresult res = kNN(*X, *Y, n, m, d, k);

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < k; j++) {
                CHECK(res.nidx[i][j] == nidx_ans[i][j]);
                CHECK(res.ndist[i][j] == ndist_ans[i][j]);
            }
        }
    }

    SUBCASE("One same point in X and Y") {
        double X[][1] = {{0.0}, {2.0}, {4.0}, {15.0}, {16.0}, {18.0}, {40.0}};
        double Y[][1] = {{0.0}, {17.0}, {30.0}};

        double ndist_ans[][3] = {{2, 4, 15}, {1, 1, 2}, {10, 12, 14}};
        int nidx_ans[][3]     = {{1, 2, 3}, {5, 4, 3}, {6, 5, 4}};

        struct knnresult res = kNN(*X, *Y, n, m, d, k);

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < k; j++) {
                CHECK(res.nidx[i][j] == nidx_ans[i][j]);
                CHECK(res.ndist[i][j] == ndist_ans[i][j]);
            }
        }
    }
}

TEST_CASE("kNN with two dimensional points") {

    int n = 2;
    int m = 3;
    int d = 2;
    int k = 2;

    double X[][2] = {{3.0, 3.0}, {5.0, 1.0}};
    double Y[][2] = {{0.0, 1.0}, {2.0, 5.0}, {4.0, 1.0}};

    double ndist_ans[][3] = {{3.60555, 5}, {2.23607, 5}, {1, 2.23607}};
    int nidx_ans[][3]     = {{0, 1}, {0, 1}, {1, 0}};

    struct knnresult res = kNN(*X, *Y, n, m, d, k);

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < k; j++) {
            CHECK(res.nidx[i][j] == nidx_ans[i][j]);
            CHECK(res.ndist[i][j] == ndist_ans[i][j]);
        }
    }
}
