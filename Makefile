# ========================================================================================
#                                  MAKEFILE MC-GPU v1.5b
# ========================================================================================

SHELL = /bin/sh
.SUFFIXES: .cu .o

# Compilers and linker:
CC = nvcc

# Program's name:
PROG = ./example_simulations/MC-GPU_v1.5b.x

# Include and library paths: 
CUDA_PATH = /usr/local/cuda/include/
CUDA_LIB_PATH = /usr/local/cuda/lib64/
CUDA_SDK_PATH = $(HOME)/cuda-samples/Common/
CUDA_SDK_LIB_PATH = $(HOME)/cuda-samples/Common/lib/x64/
OPENMPI_PATH = /usr/lib/x86_64-linux-gnu/openmpi/include/

# Base Flags (Required for both modes) 
BASE_FLAGS = -m64 -DUSING_MPI -I./ -I$(CUDA_PATH) -I$(CUDA_SDK_PATH) \
             -L$(CUDA_SDK_LIB_PATH) -L$(CUDA_LIB_PATH) \
             -lcudart -lm -lz -I$(OPENMPI_PATH) -lmpi --ptxas-options=-v

# Target-specific flags
release: CFLAGS = $(BASE_FLAGS) -O3 -use_fast_math 
debug:   CFLAGS = $(BASE_FLAGS) -g -G -O0 # -g (host) -G (device) debug symbols

# Source files: 
SRCS = MC-GPU_v1.5b.cu
RM = /bin/rm -vf 

# Default target
default: release

release: clean $(PROG)
debug: clean $(PROG)

$(PROG):
	$(CC) $(CFLAGS) $(SRCS) -o $(PROG)

clean:
	$(RM) $(PROG)