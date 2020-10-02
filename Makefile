
include make.inc

all:
	$(MAKE) -C dspl 
	$(MAKE) -C examples 
	$(MAKE) -C performance
	$(MAKE) -C verification

clean:
	$(MAKE) -C dspl clean
	$(MAKE) -C examples clean
	$(MAKE) -C performance clean
	$(MAKE) -C verification clean
	rm -f _release/*.*

clean_all:
	$(MAKE) -C dspl clean
	$(MAKE) -C dspl/blas clean
	$(MAKE) -C dspl/lapack clean
	$(MAKE) -C examples clean
	$(MAKE) -C performance clean
	$(MAKE) -C verification clean
	rm -f _release/*.*