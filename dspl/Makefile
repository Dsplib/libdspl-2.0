
include ../make.inc

# C-compiler flags
CFLAGS  = -c -fPIC -Wall -O3 -I$(INC_DIR) -Isrc -DBUILD_LIB -D$(DEF_OS)

# DSPL src and obj files list
DSPL_SRC_FILES	= $(wildcard $(DSPL_SRC_DIR)/*.c)
DSPL_OBJ_FILES	= $(addprefix $(DSPL_OBJ_DIR)/,$(notdir $(DSPL_SRC_FILES:.c=.o)))

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
$(DSPL_OBJ_DIR)/%.o:$(DSPL_SRC_DIR)/%.c
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

