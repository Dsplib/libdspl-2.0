
CC = gcc

# Define OS
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

# C-compiler flags
CFLAGS  = -c -O3 -D$(DEF_OS)

OBJFILES = test.o dspl.o

all: test.exe clean

test.exe: $(OBJFILES)
	$(CC) $(OBJFILES) -o $@ $(LFLAGS)

%.o:%.c
	$(CC) $(CFLAGS) $< -o $@ $(LFLAGS)
  
clean:
	rm -f *.o