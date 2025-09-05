#compiler settings
CXX 		  := g++
NVCC		  := nvcc
CXXFLAGS  := -std=c++17 -O2 -Iinclude
NVCCFLAGS := -std=c++17 -O2 -Iinclude

# Directories
SRC_DIR   := src
EX_DIR    := examples
BUILD_DIR := build
BIN_DIR   := bin


# Collect all sources
CPP_SRCS  := $(wildcard $(SRC_DIR)/*.cpp)
CU_SRCS   := $(wildcard $(EX_DIR)/*.cu)

# turn sources into object files
CPP_OBJS  := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(CPP_SRCS)) 
CU_OBJS   := $(patsubst $(EX_DIR)/%.cu,$(BUILD_DIR)/%.o,$(CU_SRCS)) 

# Final program
TARGET    := $(BIN_DIR)/computeAngle

#default rule
all: $(TARGET)

# Link step
$(TARGET): $(CPP_OBJS) $(CU_OBJS) | $(BIN_DIR)
	$(NVCC) -o $@ $^

# compile c++ sources
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# compile CUDA sources
$(BUILD_DIR)/%.o: $(EX_DIR)/%.cu | $(BUILD_DIR)
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

# ensure output directory exist
$(BUILD_DIR) $(BIN_DIR):
	mkdir -p $@

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)


