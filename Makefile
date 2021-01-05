MPICC=mpic++
CFLAGS= -O3 -lopenblas -std=c++17

BUILD_DIR=./build
SRC_DIR=./src
INCLUDE_DIR=./include
SOURCES := $(shell find $(SRC_DIR) -name '*.cpp')

$(info $(shell mkdir -p $(BUILD_DIR)))


default: all

mpi:
	$(MPICC) -o $(BUILD_DIR)/main -I$(INCLUDE_DIR) $(SOURCES) $(CFLAGS) 

.PHONY: clean

all: mpi test

test:
	@printf "\n** Testing\n\n"
	mpirun -np 1 ./build/main

clean:
	rm -rf $(BUILD_DIR)