MPICC=mpic++
CFLAGS= -O3 -lopenblas -std=c++17

BUILD_DIR=./build
SRC_DIR=./src
INCLUDE_DIR=./include
SOURCES := $(shell find $(SRC_DIR) -name '*.cpp')

NUM_OF_PROCESSES=8
DATA=3
VERSION=1

$(info $(shell mkdir -p $(BUILD_DIR)))


default: all

mpi:
	$(MPICC) -o $(BUILD_DIR)/main -I$(INCLUDE_DIR) $(SOURCES) $(CFLAGS) 

.PHONY: clean

all: mpi test

test:
	@printf "\n** Testing\n\n"
	mpirun -np $(NUM_OF_PROCESSES) ./build/main $(DATA) $(VERSION)

clean:
	rm -rf $(BUILD_DIR)