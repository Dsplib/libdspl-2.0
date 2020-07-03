

FORTRAN  = gfortran
OPTS     = -O3  -fPIC
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
LAPACK_DOUBLE_LIB      = bin/liblapack_double.a
LAPACK_COMPLEX_LIB      = bin/liblapack_complex.a

LAPACK_OBJ_DIR = obj

# DSPL src and obj files list
include lapack_double.inc
include lapack_complex.inc

LAPACK_DOUBLE_OBJ_FILES	  = $(addprefix $(LAPACK_OBJ_DIR)/,$(notdir $(LAPACK_DOUBLE_SRC:.f=.o)))
LAPACK_DOUBLE_OBJ_FILES1	= $(addprefix $(LAPACK_OBJ_DIR)/,$(notdir $(LAPACK_DOUBLE_SRCF:.F=.o)))
LAPACK_COMPLEX_OBJ_FILES	= $(addprefix $(LAPACK_OBJ_DIR)/,$(notdir $(LAPACK_COMPLEX_SRC:.f=.o)))
LAPACK_COMPLEX_OBJ_FILES1	= $(addprefix $(LAPACK_OBJ_DIR)/,$(notdir $(LAPACK_COMPLEX_SRCF:.F=.o)))


all: $(LAPACK_DOUBLE_LIB) $(LAPACK_COMPLEX_LIB)


$(LAPACK_DOUBLE_LIB): $(LAPACK_DOUBLE_OBJ_FILES) $(LAPACK_DOUBLE_OBJ_FILES1)
	$(ARCH) $(ARCHFLAGS) $@ $(LAPACK_DOUBLE_OBJ_FILES) $(LAPACK_DOUBLE_OBJ_FILES1)
	$(RANLIB) $@


$(LAPACK_COMPLEX_LIB): $(LAPACK_COMPLEX_OBJ_FILES) $(LAPACK_COMPLEX_OBJ_FILES1)
	$(ARCH) $(ARCHFLAGS) $@ $(LAPACK_COMPLEX_OBJ_FILES) $(LAPACK_COMPLEX_OBJ_FILES1)
	$(RANLIB) $@


$(LAPACK_OBJ_DIR)/%.o:src/%.f
	$(FORTRAN) $(OPTS) -c $< -o $@
  
$(LAPACK_OBJ_DIR)/%.o:src/%.F
	$(FORTRAN) $(OPTS) -c $< -o $@

clean:
	rm -f   $(LAPACK_OBJ_DIR)/*.o
	rm -f   $(LAPACK_DOUBLE_LIB)
	rm -f   $(LAPACK_COMPLEX_LIB)
