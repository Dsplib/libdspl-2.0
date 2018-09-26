#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

int main()
{
	void* handle;           // DSPL handle
	handle = dspl_load();   // Load DSPL function

	matrix_t r, c;

	memset(&r, 0, sizeof(matrix_t));
	memset(&c, 0, sizeof(matrix_t));

	matrix_create(&r, 4, 3, DAT_DOUBLE);
	linspace(0, 12, 12, DSPL_PERIODIC, (double*)r.dat);
	matrix_print(&r, "R", "%8.2f");

	matrix_create(&c, 2, 3, DAT_COMPLEX);
	linspace(0, 12, 12, DSPL_PERIODIC, (double*)c.dat);
	matrix_print(&c, "C", "%8.2f%+8.2fj     ");


	matrix_free(&r);
	matrix_free(&c);
	dspl_free(handle);      // free dspl handle

	// выполнить скрипт GNUPLOT для построения графиков
	// по рассчитанным данным
	return system("gnuplot gnuplot/sinc_test.plt");;
}

