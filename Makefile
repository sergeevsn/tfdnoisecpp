# Compiler [[2]]
CC = g++

# Compilation flags [[5]]
CFLAGS = -O3 -march=native -ffast-math -fopenmp -std=c++17 -I/usr/local/include -I./include

# Linker flags for FFTW3 and segyio [[6]]
LDFLAGS = -lfftw3f -lsegyio

# Build directory configuration [[4]]
BUILD_DIR = build
TARGET = $(BUILD_DIR)/tfdnoise

# Source files [[2]]
SOURCES = include/tfd.cpp tfdnoise.cpp

# Object files with build path [[4]]
OBJECTS = $(addprefix $(BUILD_DIR)/, $(SOURCES:.cpp=.o))

# Header files [[2]]
HEADERS = tfd.h

# Default target [[2]]
all: $(BUILD_DIR) $(TARGET)

# Link object files into executable [[6]]
$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -o $(TARGET) $(LDFLAGS)

# Compile source files to object files in build directory [[4]]
$(BUILD_DIR)/%.o: %.cpp $(HEADERS)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@

# Create build directory [[4]]
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Cleanup [[2]]
clean:
	rm -rf $(BUILD_DIR)/*

# Test run [[4]]
run: $(TARGET)
	./$(TARGET) params.tfd

# Phony targets [[2]]
.PHONY: all clean run
