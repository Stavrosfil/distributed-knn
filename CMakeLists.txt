cmake_minimum_required(VERSION 3.0.0)
project(knn-mpi)

# set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -lopenblas")

include(CTest)
enable_testing()

set(EXEC main)

include_directories(include)

file(GLOB SOURCES src/*.cpp)
add_executable(${EXEC} ${SOURCES})

# add_executable(knn-mpi src/main.cpp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

