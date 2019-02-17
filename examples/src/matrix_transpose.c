#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

int main()
{
	void* handle;           // DSPL handle
	handle = dspl_load();   // Load DSPL function

	matrix_t a, b, c;

	memset(&a, 0, sizeof(matrix_t));
	memset(&b, 0, sizeof(matrix_t));
	memset(&c, 0, sizeof(matrix_t));

	matrix_create(&a, 4, 3, DAT_DOUBLE);
	linspace(0, 12, 12, DSPL_PERIODIC, (double*)a.dat);
	matrix_print(&a, "A", "%8.2f");

	matrix_transpose(&a, &b);
	matrix_print(&b, "B=A^T", "%8.2f");

	matrix_transpose_hermite(&a, &b);
	matrix_print(&b, "B=A^H", "%8.2f");



	matrix_create(&c, 2, 3, DAT_COMPLEX);
	linspace(0, 12, 12, DSPL_PERIODIC, (double*)c.dat);
	matrix_print(&c, "C", "%8.2f%+8.2fj     ");

	matrix_transpose(&c, &b);
	matrix_print(&b, "B=C^T", "%8.2f%+8.2fj     ");


	matrix_transpose_hermite(&c, &b);
	matrix_print(&b, "B=C^H", "%8.2f%+8.2fj     ");


	matrix_free(&a);
	matrix_free(&b);
	matrix_free(&c);
	dspl_free(handle);      // free dspl handle

	// выполнить скрипт GNUPLOT для построения графиков
	// по рассчитанным данным
	return system("gnuplot gnuplot/sinc_test.plt");;
}

