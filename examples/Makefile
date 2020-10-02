
include ../make.inc


SRC_DIR = src
OBJ_DIR = obj

# C-compiler flags
CFLAGS  = -c -O3 -I$(INC_DIR) -D$(DEF_OS)


DSPL_C = ../include/dspl.c
DSPL_O = obj/dspl.o


# DSPL src and obj files list
SRC_FILES	= $(wildcard $(SRC_DIR)/*.c)
OBJ_FILES	= $(addprefix $(OBJ_DIR)/,$(notdir $(SRC_FILES:.c=.o)))
EXE_FILES	= $(addprefix $(EXAMPLE_BIN_DIR)/,$(notdir $(SRC_FILES:.c=.exe)))


all: $(EXE_FILES)

$(EXAMPLE_BIN_DIR)/%.exe:$(OBJ_DIR)/%.o $(DSPL_O)
	$(CC)  $< $(DSPL_O) -o $@ $(LFLAGS)

$(OBJ_DIR)/%.o:$(SRC_DIR)/%.c
	$(CC) $(CFLAGS) $< -o $@ $(LFLAGS)

$(DSPL_O):$(DSPL_C)
	$(CC) $(CFLAGS) $(DSPL_C) -o $(DSPL_O)

clean:
	rm -f   $(OBJ_DIR)/*.o
	rm -f   $(EXAMPLE_BIN_DIR)/*.exe


