MPICC=mpic++
CFLAGS= -O3 -w -lopenblas

BUILD_DIR=build
SRC_DIR=src
INCLUDE_DIR=include
SOURCES := $(shell find $(SRC_DIR) -name '*.cpp')

$(info $(shell mkdir -p $(BUILD_DIR)))


default: all

mpi:
	$(MPICC) $(CFLAGS) -o $(BUILD_DIR)/main -I$(INCLUDE_DIR) $(SOURCES)

.PHONY: clean

all: mpi test

test:
	@printf "\n** Testing\n\n"
	./build/main

clean:
	rm -rf $(BUILD_DIR)