
include ../make.inc

# C-compiler flags
CFLAGS  = -c -fPIC -Wall -O3 -I$(INC_DIR) -DBUILD_LIB -D$(DEF_OS)

# DSPL src and obj files list
DSPL_SRC_FILES	= $(wildcard $(DSPL_SRC_DIR)/*.c)
DSPL_OBJ_FILES	= $(addprefix $(DSPL_OBJ_DIR)/,$(notdir $(DSPL_SRC_FILES:.c=.o)))

all: $(RELEASE_DIR)/$(LIB_NAME)



$(RELEASE_DIR)/$(LIB_NAME): $(DSPL_OBJ_FILES)  $(BLAS_LIB_NAME) $(LAPACK_DOUBLE_LIB_NAME) $(LAPACK_COMPLEX_LIB_NAME)
	$(CC) -shared -o $(RELEASE_DIR)/$(LIB_NAME)  $(DSPL_OBJ_FILES) -lm  -L$(LAPACK_RELEASE_DIR)  -llapack_complex -llapack_double -L$(BLAS_RELEASE_DIR) -lblas -lgfortran -lquadmath


$(DSPL_OBJ_DIR)/%.o:$(DSPL_SRC_DIR)/%.c
	$(CC) $(CFLAGS)  $< -o $@ -lm


clean:
	rm -f   $(DSPL_OBJ_DIR)/*.o
	rm -f   $(RELEASE_DIR)/$(LIB_NAME)

