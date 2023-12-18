
include ../make.inc

# C-compiler flags
CFLAGS  = -c -fPIC -Wall -O3 -I$(INC_DIR) -Isrc -DBUILD_LIB -D$(DEF_OS)

# DSPL src list according to pathes
ARRAY_C         = $(wildcard $(DSPL_SRC_DIR)/array/*.c)
CONV_C          = $(wildcard $(DSPL_SRC_DIR)/convolution/*.c)
DFT_C           = $(wildcard $(DSPL_SRC_DIR)/dft/*.c)
FILTER_C        = $(wildcard $(DSPL_SRC_DIR)/filter_design/*.c)
GNUPLOT_C       = $(wildcard $(DSPL_SRC_DIR)/gnuplot/*.c)
INOUT_C         = $(wildcard $(DSPL_SRC_DIR)/inout/*.c)
MATH_ELLIPJ_C   = $(wildcard $(DSPL_SRC_DIR)/math_ellipj/*.c)
MATH_MATRIX_C   = $(wildcard $(DSPL_SRC_DIR)/math_matrix/*.c)
MATH_POLY_C     = $(wildcard $(DSPL_SRC_DIR)/math_poly/*.c)
MATH_STAT_C     = $(wildcard $(DSPL_SRC_DIR)/math_stat/*.c)
MATH_STD_C      = $(wildcard $(DSPL_SRC_DIR)/math_std/*.c)
PSD_C           = $(wildcard $(DSPL_SRC_DIR)/psd/*.c)
RANDOM_C        = $(wildcard $(DSPL_SRC_DIR)/randomgen/*.c)
RESAMPLING_C    = $(wildcard $(DSPL_SRC_DIR)/resampling/*.c)
SIGNALS_C       = $(wildcard $(DSPL_SRC_DIR)/signals/*.c)
TYPES_C         = $(wildcard $(DSPL_SRC_DIR)/types/*.c)

# DSPL OBJ list according to pathes
ARRAY_O         = $(addprefix $(DSPL_OBJ_DIR)/,$(notdir $(ARRAY_C:.c=.o)))
CONV_O          = $(addprefix $(DSPL_OBJ_DIR)/,$(notdir $(CONV_C:.c=.o)))
DFT_O           = $(addprefix $(DSPL_OBJ_DIR)/,$(notdir $(DFT_C:.c=.o)))
FILTER_O        = $(addprefix $(DSPL_OBJ_DIR)/,$(notdir $(FILTER_C:.c=.o)))
GNUPLOT_O       = $(addprefix $(DSPL_OBJ_DIR)/,$(notdir $(GNUPLOT_C:.c=.o)))
INOUT_O         = $(addprefix $(DSPL_OBJ_DIR)/,$(notdir $(INOUT_C:.c=.o)))
MATH_ELLIPJ_O   = $(addprefix $(DSPL_OBJ_DIR)/,$(notdir $(MATH_ELLIPJ_C:.c=.o)))
MATH_MATRIX_O   = $(addprefix $(DSPL_OBJ_DIR)/,$(notdir $(MATH_MATRIX_C:.c=.o)))
MATH_POLY_O     = $(addprefix $(DSPL_OBJ_DIR)/,$(notdir $(MATH_POLY_C:.c=.o)))
MATH_STAT_O     = $(addprefix $(DSPL_OBJ_DIR)/,$(notdir $(MATH_STAT_C:.c=.o)))
MATH_STD_O      = $(addprefix $(DSPL_OBJ_DIR)/,$(notdir $(MATH_STD_C:.c=.o)))
PSD_O           = $(addprefix $(DSPL_OBJ_DIR)/,$(notdir $(PSD_C:.c=.o)))
RANDOM_O        = $(addprefix $(DSPL_OBJ_DIR)/,$(notdir $(RANDOM_C:.c=.o)))
RESAMPLING_O    = $(addprefix $(DSPL_OBJ_DIR)/,$(notdir $(RESAMPLING_C:.c=.o)))
SIGNALS_O       = $(addprefix $(DSPL_OBJ_DIR)/,$(notdir $(SIGNALS_C:.c=.o)))
TYPES_O         = $(addprefix $(DSPL_OBJ_DIR)/,$(notdir $(TYPES_C:.c=.o)))


DSPL_OBJ_FILES  = $(ARRAY_O)\
                  $(CONV_O)\
                  $(DFT_O)\
                  $(FILTER_O)\
                  $(GNUPLOT_O)\
                  $(INOUT_O)\
                  $(MATH_ELLIPJ_O)\
                  $(MATH_MATRIX_O)\
                  $(MATH_POLY_O)\
                  $(MATH_STAT_O)\
                  $(MATH_STD_O)\
                  $(PSD_O)\
                  $(RANDOM_O)\
                  $(RESAMPLING_O)\
                  $(SIGNALS_O)\
                  $(TYPES_O)\


all: $(RELEASE_DIR)/$(LIB_NAME)\
     $(EXAMPLE_BIN_DIR)/$(LIB_NAME)\
     $(PERFORMANCE_BIN_DIR)/$(LIB_NAME)\
     $(VERIFICATION_BIN_DIR)/$(LIB_NAME)\
     $(RELEASE_DIR)/dspl.c\
     $(RELEASE_DIR)/dspl.h

#Build libdspl.dll or libdspl.so
$(RELEASE_DIR)/$(LIB_NAME): $(DSPL_OBJ_FILES)  $(BLAS_LIB_NAME) $(LAPACK_DOUBLE_LIB_NAME) $(LAPACK_COMPLEX_LIB_NAME)
	$(CC) -shared -o $(RELEASE_DIR)/$(LIB_NAME)  $(DSPL_OBJ_FILES) -lm  -L$(LAPACK_RELEASE_DIR) -llapack_complex -llapack_double -L$(BLAS_RELEASE_DIR) -lblas -lgfortran -lquadmath



#Compile libdspl obj files from c sources
$(DSPL_OBJ_DIR)/%.o:$(DSPL_SRC_DIR)/array/%.c
	$(CC) $(CFLAGS)  $< -o $@ -lm 
$(DSPL_OBJ_DIR)/%.o:$(DSPL_SRC_DIR)/convolution/%.c
	$(CC) $(CFLAGS)  $< -o $@ -lm 
$(DSPL_OBJ_DIR)/%.o:$(DSPL_SRC_DIR)/dft/%.c
	$(CC) $(CFLAGS)  $< -o $@ -lm 
$(DSPL_OBJ_DIR)/%.o:$(DSPL_SRC_DIR)/filter_design/%.c
	$(CC) $(CFLAGS)  $< -o $@ -lm 
$(DSPL_OBJ_DIR)/%.o:$(DSPL_SRC_DIR)/gnuplot/%.c
	$(CC) $(CFLAGS)  $< -o $@ -lm 
$(DSPL_OBJ_DIR)/%.o:$(DSPL_SRC_DIR)/inout/%.c
	$(CC) $(CFLAGS)  $< -o $@ -lm 
$(DSPL_OBJ_DIR)/%.o:$(DSPL_SRC_DIR)/math_ellipj/%.c
	$(CC) $(CFLAGS)  $< -o $@ -lm 
$(DSPL_OBJ_DIR)/%.o:$(DSPL_SRC_DIR)/math_matrix/%.c
	$(CC) $(CFLAGS)  $< -o $@ -lm 
$(DSPL_OBJ_DIR)/%.o:$(DSPL_SRC_DIR)/math_poly/%.c
	$(CC) $(CFLAGS)  $< -o $@ -lm 
$(DSPL_OBJ_DIR)/%.o:$(DSPL_SRC_DIR)/math_stat/%.c
	$(CC) $(CFLAGS)  $< -o $@ -lm 
$(DSPL_OBJ_DIR)/%.o:$(DSPL_SRC_DIR)/math_std/%.c
	$(CC) $(CFLAGS)  $< -o $@ -lm 
$(DSPL_OBJ_DIR)/%.o:$(DSPL_SRC_DIR)/psd/%.c
	$(CC) $(CFLAGS)  $< -o $@ -lm 
$(DSPL_OBJ_DIR)/%.o:$(DSPL_SRC_DIR)/randomgen/%.c
	$(CC) $(CFLAGS)  $< -o $@ -lm 
$(DSPL_OBJ_DIR)/%.o:$(DSPL_SRC_DIR)/resampling/%.c
	$(CC) $(CFLAGS)  $< -o $@ -lm 
$(DSPL_OBJ_DIR)/%.o:$(DSPL_SRC_DIR)/signals/%.c
	$(CC) $(CFLAGS)  $< -o $@ -lm 
$(DSPL_OBJ_DIR)/%.o:$(DSPL_SRC_DIR)/types/%.c
	$(CC) $(CFLAGS)  $< -o $@ -lm 


#Copy libdspl.dll to the examples "bin" folder 
$(EXAMPLE_BIN_DIR)/$(LIB_NAME):$(RELEASE_DIR)/$(LIB_NAME)
	cp $(RELEASE_DIR)/$(LIB_NAME) $(EXAMPLE_BIN_DIR)/$(LIB_NAME)

#Copy libdspl.dll to the performance testing "bin" folder
$(PERFORMANCE_BIN_DIR)/$(LIB_NAME):$(RELEASE_DIR)/$(LIB_NAME)
	cp $(RELEASE_DIR)/$(LIB_NAME) $(PERFORMANCE_BIN_DIR)/$(LIB_NAME)

#Copy libdspl.dll to the verification "bin" folder
$(VERIFICATION_BIN_DIR)/$(LIB_NAME):$(RELEASE_DIR)/$(LIB_NAME)
	cp $(RELEASE_DIR)/$(LIB_NAME) $(VERIFICATION_BIN_DIR)/$(LIB_NAME)

#make BLAS library
$(BLAS_LIB_NAME):
	$(MAKE) -C blas

#make LAPACK library for the real data
$(LAPACK_DOUBLE_LIB_NAME):
	$(MAKE) -C lapack

#make LAPACK library for the complex data
$(LAPACK_COMPLEX_LIB_NAME):
	$(MAKE) -C lapack

#Copy dspl.c source to the release folder 
$(RELEASE_DIR)/dspl.c:
	cp $(INC_DIR)/dspl.c $(RELEASE_DIR)/dspl.c

#Copy dspl.h header to the release folder
$(RELEASE_DIR)/dspl.h:
	cp $(INC_DIR)/dspl.h $(RELEASE_DIR)/dspl.h


clean:
	rm -f   $(DSPL_OBJ_DIR)/*.o
	rm -f   $(RELEASE_DIR)/*.a
	rm -f   $(RELEASE_DIR)/*.def

