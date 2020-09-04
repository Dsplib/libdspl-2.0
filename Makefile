
include make.inc

all:
	$(MAKE) -C blas
	$(MAKE) -C lapack
	$(MAKE) -C dspl 
	$(MAKE) -C examples 
	$(MAKE) -C performance
	$(MAKE) -C verification

clean:
	$(MAKE) -C dspl clean
	$(MAKE) -C examples clean
	$(MAKE) -C performance clean
	$(MAKE) -C verification clean