# Distributed K-Nearest-Neighbors using MPI and a Vantage Point Tree

## Datasets used for benchmarking

- [corel](https://archive.ics.uci.edu/ml/datasets/Corel+Image+Features)
- [fma](https://archive.ics.uci.edu/ml/datasets/FMA%3A+A+Dataset+For+Music+Analysis)
- [miniboone](https://archive.ics.uci.edu/ml/datasets/MiniBooNE+particle+identification)
- [tv_news_com](https://archive.ics.uci.edu/ml/datasets/TV+News+Channel+Commercial+Detection+Dataset])

## Build instructions

You will need an mpi compatible compiler (`mpic++`) to compile and run (`mpirun`) the application.

You can find more details about compilation inside the `Makefile`.

```bash
make
```

## Instructions for dataset paths

Place each file named `$fileName` that belongs to the `$dataset` in the directory with the following structure:

```
$repoRootDirectory/dataset/$dataset/$fileName
```

For example, for the file `ColorHistogram.asc` from `corel` dataset use

```
./dataset/corel/ColorHistogram.asc
```

Note: Some files are cleaned using a python script before processing, beware.