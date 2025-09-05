NVCC=nvcc
CFLAGS = -O3 ## optimize make


#Directories
SRC_DIR = src
BIN_DIR = bin


#programs to build
#PROGRAMS = conductivity rate
## automatically detect all subdirectories`
PROGRAMS := $(notdir $(wildcard $(SRC_DIR)/*))

define SRC_FILES
$(wildcard $(SRC_DIR)/$1/*.cu)
endef

## create map: executable depends on its source files
define PROGRAM_RULE
$(BIN_DIR)/$1: $$(call SRC_FILES,$1)
	@echo "Building $1 ..."
	@mkdir -p $(BIN_DIR)
	$$(NVCC) $$(CFLAGS) $$(call SRC_FILES,$1) -o $$@
endef

# Target executable
TARGETS = $(patsubst %, $(BIN_DIR)/%, $(PROGRAMS))

.PHONY: all clean

all: $(TARGETS)

# one build rule per program
$(foreach prog,$(PROGRAMS), $(eval $(call PROGRAM_RULE,$(prog))))

# allow "make program1" instead of "make bin/program1"
$(PROGRAMS): %: $(BIN_DIR)/%
	@echo "Done building $@"

clean:
	rm -rf $(BIN_DIR)


