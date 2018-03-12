



ifeq ($(OS),Windows_NT)
    MAKE = mingw32-make  
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
        MAKE = make       
    endif
endif

include Makefile.dirs




all:
	$(MAKE) -f Makefile.dspl 
	cp  -r include/dspl.h release/include/dspl.h
	cp  -r include/dspl.c release/include/dspl.c
	$(MAKE) -f Makefile.test
	cp  -r $(RELEASE_DIR)/$(DSPL_LIBNAME) test/bin/$(DSPL_LIBNAME)
 
	
	

clean:
	$(MAKE) -f Makefile.dspl  clean
	$(MAKE) -f Makefile.test  clean
	




	

