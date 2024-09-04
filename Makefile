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
PARALLEL_TESTS_SRC := $(SRC_DIR)/parallel_tests.cpp
SEQUENTIAL_TESTS_SRC := $(SRC_DIR)/sequential_tests.cpp

# Executables
SEQUENTIAL_EXEC := $(BIN_DIR)/sequential
PARALLEL_EXEC := $(BIN_DIR)/parallel
PARALLEL_TESTS_EXEC := $(BIN_DIR)/parallel_tests
SEQUENTIAL_TESTS_EXEC := $(BIN_DIR)/sequential_tests

# Default target
all: $(SEQUENTIAL_EXEC) $(PARALLEL_EXEC) $(PARALLEL_TESTS_EXEC) $(SEQUENTIAL_TESTS_EXEC)

# Targets for individual executables
sequential: $(SEQUENTIAL_EXEC)

parallel: $(PARALLEL_EXEC)

parallel_tests: $(PARALLEL_TESTS_EXEC)

sequential_tests: $(SEQUENTIAL_TESTS_EXEC)

$(SEQUENTIAL_EXEC): $(SEQUENTIAL_SRC) | $(BIN_DIR)
	$(CXX) $(SEQUENTIAL_FLAGS) $< -o $@ $(SEQUENTIAL_LDFLAGS)

$(PARALLEL_EXEC): $(PARALLEL_SRC) | $(BIN_DIR)
	$(CXX) $(PARALLEL_FLAGS) $< -o $@ $(PARALLEL_LDFLAGS)

$(PARALLEL_TESTS_EXEC): $(PARALLEL_TESTS_SRC) | $(BIN_DIR)
	$(CXX) $(PARALLEL_FLAGS) $< -o $@ $(PARALLEL_LDFLAGS)

$(SEQUENTIAL_TESTS_EXEC): $(SEQUENTIAL_TESTS_SRC) | $(BIN_DIR)
	$(CXX) $(SEQUENTIAL_FLAGS) $< -o $@ $(SEQUENTIAL_LDFLAGS)

$(BIN_DIR):
	mkdir -p $@

clean:
	rm -rf $(BIN_DIR)

.PHONY: all sequential parallel parallel_tests sequential_tests clean
