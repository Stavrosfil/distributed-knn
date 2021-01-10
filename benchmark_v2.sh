#!/bin/bash

#SBATCH --job-name=pds_project
# This is SBATCH --partition=batch

# Number of MPI tasks
#SBATCH -n 40
#
# Number of tasks per node
#SBATCH --tasks-per-node=20
#
# Runtime of this jobs is less then 12 hours.
#SBATCH --time=00:01:00

module load gcc
module load openmpi
module load openblas

processors_nums=(1 2 4 8)
datasets=(0 1 2 3)

mpic++ -o ./build/main -I./include ./src/main.cpp -O3 -lopenblas -std=c++17

for i in ${datasets[@]}; do
    for j in ${processors_nums[@]}; do
        printf "\nVERSION: v2 | PROCESSORS: $j\n"
        srun -n $j ./build/main $i 2
        echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    done
done
