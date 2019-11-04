
include make.inc

all:
	$(MAKE) -C blas
	$(MAKE) -C lapack
	$(MAKE) -C dspl 
	$(MAKE) -C examples 


clean:
	$(MAKE) -C dspl clean
	$(MAKE) -C examples clean