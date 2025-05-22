# Makefile for compiling main.c

# Compiler to use
CC = gcc

# Compiler flags
CFLAGS = -Wall -Wextra -O3 -ffast-math

# The target directory for output
BUILD_DIR = build

# The target executable
TARGET = $(BUILD_DIR)/ntt_bench

# Source files
SRC = main.c

# Create build directory if it doesn't exist
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Build target
all: $(BUILD_DIR) $(TARGET)

# Linking the object file to create the executable
$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC)

# Clean up the compiled files
clean:
	rm -rf $(BUILD_DIR)

.PHONY: all clean