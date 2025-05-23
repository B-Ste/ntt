# Define the compiler to use
CC = gcc

# Define compiler flags
CFLAGS = -Wall -Wextra -O3 -ffast-math

# Define the source file
SRC = main.c

# Define the build directory and the target executable name
BUILD_DIR = build
TARGET = $(BUILD_DIR)/ntt_bench

# Targets to create
.PHONY: all with_bar without_bar clean

all: without_bar with_bar

# Target for compiling without the -D BAR flag
without_bar: $(TARGET)

$(TARGET): $(SRC)
	mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) -o $@ $^

# Target for compiling with the -D BAR flag
with_bar: $(TARGET)_bar

$(TARGET)_bar: $(SRC)
	mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) -D BAR -o $@ $^

# Clean up the build directory
clean:
	rm -rf $(BUILD_DIR)
