#include <stdio.h>
#include <string.h>
#include <math.h>
#include "dspl.h"

#define N  1000
#define T  4
#define TAU 2
#define K 41

int main(int argc, char* argv[])
{
	double t[N];	// время
	double s[N];	// исходный сигнал
	double w[K];	// массив частоты
	complex_t S[K]; // спектр
	double mag[K];	// амплитудный спектр
	double phi[K];	// фазовый спектр
	void* handle;
	int k, n;

	handle = dspl_load();
	if(!handle)
	{
		printf("cannot to load libdspl!\n");
		return 0;
	}

	// заполняем массив времени для одного периода сигнала
	linspace(-T*0.5, T*0.5, N, DSPL_SYMMETRIC, t);

	// сигнал
	signal_pimp(t, N, 2.0, TAU, 0.0, T, s);

	// расчет комплексного спектра
	fourier_series_dec(t, s, N, T, K, w, S);

	// Амплитудный и фазовый спектры
	for(n = 0; n < K; n++)
	{
		mag[n] = ABS(S[n]);
		phi[n] = atan2(IM(S[n]), RE(S[n]));
	}

	// Сохранение амплитудного спектра в файл
	writetxt(w, mag, K, "dat/fourier_series_pimp_mag.txt");

	// Сохранение фазового спектра в файл
	writetxt(w, phi, K, "dat/fourier_series_pimp_phi.txt");


	// remember to free the resource
	dspl_free(handle);

	return system("gnuplot -p  gnuplot/fourier_series_pimp.plt");
}
