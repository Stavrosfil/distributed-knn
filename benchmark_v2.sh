#!/bin/bash

processors_nums=(1 2 4 8)
datasets=(0 1 2 3)

mpic++ -o ./build/main -I./include ./src/main.cpp -O3 -lopenblas -std=c++17

for i in ${datasets[@]}
do
    for j in ${processors_nums[@]}
    do
        printf "\nVERSION: v2 | PROCESSORS: $j\n"
        mpirun -np $j ./build/main $i 2
        echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    done
done
