#ifndef __UTILS_H__
#define __UTILS_H__

#include <iostream>
#include <stdlib.h>

namespace prt {

void rowMajor(double* a, int n, int m) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++)
            std::cout << a[i * m + j] << " ";
        std::cout << std::endl;
    }
}

void twoDim(int* a, int n, int m) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++)
            std::cout << a[i * m + j] << " ";
        std::cout << std::endl;
    }
}

void twoDim(double* a, int n, int m) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++)
            std::cout << a[i * m + j] << " ";
        std::cout << std::endl;
    }
}

} // namespace prt
#endif // __UTILS_H__