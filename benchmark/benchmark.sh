#!/bin/bash

#SBATCH --job-name=pds_project
# This is SBATCH --partition=batch

# Number of MPI tasks
#SBATCH -n 60
#
# Number of tasks per node
#SBATCH --tasks-per-node=20
#
# Runtime of this jobs is less then 12 hours.
#SBATCH --time=00:10:00

module load gcc
module load openmpi
module load openblas

processors_nums=(1 2 4 8 16 32 48 60)
datasets=(4)
version=(2)

mpic++ -o ./build/main -I./include ./src/main.cpp -O3 -lopenblas -std=c++17

for i in ${datasets[@]}; do
    for j in ${processors_nums[@]}; do
        printf "\nVERSION: v$version | PROCESSORS: $j\n"
        srun -n $j ./build/main $i $version
        echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    done
done
