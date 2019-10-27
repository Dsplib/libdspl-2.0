
include make.inc

all:
	$(MAKE) -C blas
	$(MAKE) -C lapack
	$(MAKE) -C dspl 



clean:
	$(MAKE) -C dspl clean
	$(MAKE) -C blas clean
	$(MAKE) -C lapack clean
