
include make.inc

all:
	$(MAKE) -f Makefile.dspl 
	cp  -r include/dspl.h release/include/dspl.h
	cp  -r include/dspl.c release/include/dspl.c
	$(MAKE) -f Makefile.verif
	cp  -r $(RELEASE_DIR)/$(DSPL_LIBNAME) verif/bin/$(DSPL_LIBNAME)
	$(MAKE) -f Makefile.examples
	cp  -r $(RELEASE_DIR)/$(DSPL_LIBNAME) examples/bin/$(DSPL_LIBNAME)

	
	

clean:
	$(MAKE) -f Makefile.dspl  clean
	$(MAKE) -f Makefile.verif  clean
	$(MAKE) -f Makefile.examples  clean
	rm -f $(BLAS_SRC_DIR)/*.o
	rm -f $(LAPACK_SRC_DIR)/*.o


clean_all:
	$(MAKE) clean
	rm -f $(BLAS_LIB_DIR)/*.a
	rm -f $(LAPACK_LIB_DIR)/*.a

