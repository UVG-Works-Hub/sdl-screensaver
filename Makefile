CXX := g++
CXXFLAGS := -Wall -Wextra -g -std=c++11
LDFLAGS := -lSDL2 -lm -fopenmp

# Directories
SRC_DIR := src
BIN_DIR := bin

# Source files
SEQUENTIAL_SRC := $(SRC_DIR)/sequential.cpp
PARALLEL_SRC := $(SRC_DIR)/parallel.cpp

# Executables
SEQUENTIAL_EXEC := $(BIN_DIR)/mandelbrot_sequential
PARALLEL_EXEC := $(BIN_DIR)/mandelbrot_parallel

# Default target
all: $(SEQUENTIAL_EXEC) $(PARALLEL_EXEC)

# Targets for individual executables
sequential: $(SEQUENTIAL_EXEC)

parallel: $(PARALLEL_EXEC)

$(SEQUENTIAL_EXEC): $(SEQUENTIAL_SRC) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

$(PARALLEL_EXEC): $(PARALLEL_SRC) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

$(BIN_DIR):
	mkdir -p $@

clean:
	rm -rf $(BIN_DIR)

.PHONY: all sequential parallel clean
