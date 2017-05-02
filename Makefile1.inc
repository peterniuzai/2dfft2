# The definitions in this file should be changed to match your system

CUDA_PATH ?= /usr/local/cuda-8.0
INSTALL_DIR = $(HOME)/linux_64
#INSTALL_DIR = /home/nch/opt
CUDA_DIR   ?= $(CUDA_PATH)
#THRUST_DIR ?= $(CUDA_DIR)/include
#THRUST_DIR  ?= /usr/local/cuda-8.0/include/thrust
THRUST_DIR  ?= $(HOME)/opt/thrust-1.6

# LIB_ARCH   = lib
LIB_ARCH   = lib64
# GPU_ARCH   = sm_13
# GPU_ARCH   = sm_20
GPU_ARCH   = sm_61

GCC        = gcc
GXX        = g++
AR         = ar
NVCC       = $(CUDA_DIR)/bin/nvcc #-Xptxas -abi=no
DOXYGEN    = doxygen
RM         = rm
ECHO       = echo
MKFLAGS    = 

DEDISP_DEBUG = 0

#.SILENT :
