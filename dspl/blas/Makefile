

FORTRAN  = g++
OPTS     = -O3 -fPIC
DRVOPTS  = $(OPTS)
NOOPT    =
LOADOPTS =
#
#  The archiver and the flag(s) to use when building archive (library)
#  If you system has no ranlib, set RANLIB = echo.
#
ARCH     = ar
ARCHFLAGS= cr
RANLIB   = ranlib
#
#  The location and name of the Reference BLAS library.
#
BLASLIB      = bin/libblas.a


BLAS_OBJ_DIR = obj

# DSPL src and obj files list
BLAS_SRC_FILES	= $(wildcard src/*.f)
BLAS_OBJ_FILES	= $(addprefix $(BLAS_OBJ_DIR)/,$(notdir $(BLAS_SRC_FILES:.f=.o)))

all: $(BLASLIB)


$(BLASLIB): $(BLAS_OBJ_FILES)
	$(ARCH) $(ARCHFLAGS) $@ $(BLAS_OBJ_FILES)
	$(RANLIB) $@


$(BLAS_OBJ_DIR)/%.o:src/%.f
	$(FORTRAN) $(OPTS) -c $< -o $@


clean:
	rm -f   $(BLAS_OBJ_DIR)/*.o
	rm -f   $(BLASLIB)

