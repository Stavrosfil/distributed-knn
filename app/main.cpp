#ifdef ENABLE_DOCTEST_IN_LIBRARY
#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"
#endif

#include <iostream>
#include <stdlib.h>

#include "exampleConfig.h"
#include "example.h"

/*
 * Simple main program that demontrates how access
 * CMake definitions (here the version number) from source code.
 */
int main() {
    std::cout << "This is a test to check that everything is working as intended." << std::endl;

    // Bring in the dummy class from the example source,
    // just to show that it is accessible from main.cpp.
    Dummy d = Dummy();
    return d.doSomething() ? 0 : -1;
}
