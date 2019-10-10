CC       = gcc
FORTRAN  = gfortran
AR       = ar


LIB_DIR    = dspl

# DSPL source and obj file path
DSPL_SRC_DIR = $(LIB_DIR)/dspl_src
DSPL_OBJ_DIR = $(LIB_DIR)/dspl_obj

# BLAS source path
BLAS_LIB_DIR = $(LIB_DIR)/libblas/
BLAS_SRC_DIR = $(BLAS_LIB_DIR)/SRC
BLAS_LIB_NAME = $(BLAS_LIB_DIR)/libblas.a


# LAPACK source path
LAPACK_LIB_DIR = $(LIB_DIR)/liblapack/
LAPACK_SRC_DIR = $(LAPACK_LIB_DIR)/SRC
LAPACK_LIB_NAME = $(LAPACK_LIB_DIR)/liblapack.a

#common source for DSPL and examples
COMMON_SRC_DIR = $(COMMON_DIR)/src

INC_DIR = include
RELEASE_DIR = release/lib

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