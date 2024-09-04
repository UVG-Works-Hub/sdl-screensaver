CXX := g++
COMMON_FLAGS := -Wall -Wextra -std=c++17
OPTIMIZATION_FLAGS := -O3 -march=native -mtune=native -ffast-math -funroll-loops
DEBUG_FLAGS := -g
OPENMP_FLAGS := -fopenmp

SEQUENTIAL_FLAGS := $(COMMON_FLAGS) $(OPTIMIZATION_FLAGS) $(DEBUG_FLAGS)
PARALLEL_FLAGS := $(COMMON_FLAGS) $(OPTIMIZATION_FLAGS) $(DEBUG_FLAGS) $(OPENMP_FLAGS) -mavx2

SEQUENTIAL_LDFLAGS := -lSDL2 -lm
PARALLEL_LDFLAGS := -lSDL2 -lm $(OPENMP_FLAGS)

# Directories
SRC_DIR := src
BIN_DIR := bin

# Source files
SEQUENTIAL_SRC := $(SRC_DIR)/sequential.cpp
PARALLEL_SRC := $(SRC_DIR)/parallel.cpp

# Executables
SEQUENTIAL_EXEC := $(BIN_DIR)/sequential
PARALLEL_EXEC := $(BIN_DIR)/parallel

# Default target
all: $(SEQUENTIAL_EXEC) $(PARALLEL_EXEC)

# Targets for individual executables
sequential: $(SEQUENTIAL_EXEC)

parallel: $(PARALLEL_EXEC)

$(SEQUENTIAL_EXEC): $(SEQUENTIAL_SRC) | $(BIN_DIR)
	$(CXX) $(SEQUENTIAL_FLAGS) $< -o $@ $(SEQUENTIAL_LDFLAGS)

$(PARALLEL_EXEC): $(PARALLEL_SRC) | $(BIN_DIR)
	$(CXX) $(PARALLEL_FLAGS) $< -o $@ $(PARALLEL_LDFLAGS)

$(BIN_DIR):
	mkdir -p $@

clean:
	rm -rf $(BIN_DIR)

.PHONY: all sequential parallel clean
