CC       = gcc
FORTRAN  = g++
AR       = ar

# DSPL source and obj file path
DSPL_SRC_DIR = src
DSPL_OBJ_DIR = obj

# BLAS source path
BLAS_LIB_DIR = blas
BLAS_SRC_DIR = $(BLAS_LIB_DIR)/src
BLAS_RELEASE_DIR = $(BLAS_LIB_DIR)/bin
BLAS_LIB_NAME = $(BLAS_RELEASE_DIR)/libblas.a

# LAPACK source path
LAPACK_LIB_DIR  =   lapack
LAPACK_RELEASE_DIR = $(LAPACK_LIB_DIR)/bin

LAPACK_DOUBLE_SRC_DIR  = $(LAPACK_LIB_DIR)/src
LAPACK_DOUBLE_LIB_NAME = $(LAPACK_RELEASE_DIR)/liblapack_double.a

LAPACK_COMPLEX_SRC_DIR  = $(LAPACK_LIB_DIR)/src
LAPACK_COMPLEX_LIB_NAME = $(LAPACK_RELEASE_DIR)/liblapack_complex.a



INC_DIR = ../include
RELEASE_DIR = ../_release
EXAMPLE_BIN_DIR = ../examples/bin
VERIFICATION_BIN_DIR = ../verification/bin
PERFORMANCE_BIN_DIR = ../performance/bin

ifeq ($(OS),Windows_NT)
	DSPL_LIBNAME = libdspl.dll
	DEF_OS = WIN_OS
	LFLAGS = -lm
else
	UNAME_S := $(shell uname -s)
	UNAME_P := $(shell uname -p)
	ifeq ($(UNAME_S),Linux)
		DSPL_LIBNAME = libdspl.so
		DEF_OS = LINUX_OS
		LFLAGS = -lm -ldl
	else ifeq ($(UNAME_S),Darwin)
		DSPL_LIBNAME = libdspl.so
		DEF_OS = LINUX_OS
		LFLAGS = -lm -ldl
	endif
endif

LIB_NAME = $(DSPL_LIBNAME)


ifeq ($(OS),Windows_NT)
    MAKE = mingw32-make  
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
        MAKE = make
    else ifeq ($(UNAME_S),Darwin)
    	MAKE = make
    endif
endif