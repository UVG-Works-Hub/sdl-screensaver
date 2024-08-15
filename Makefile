CXX := g++
CXXFLAGS := -Wall -Wextra -g -std=c++11
LDFLAGS := -lSDL2 -lm -fopenmp

# Directories
SRC_DIR := src
BIN_DIR := bin

# Source and executable
SRC := $(SRC_DIR)/main.cpp
EXEC := $(BIN_DIR)/screensaver

# Targets
all: $(EXEC)

$(EXEC): $(SRC) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

$(BIN_DIR):
	mkdir -p $@

clean:
	rm -rf $(BIN_DIR)

.PHONY: all clean