#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"


#define N  200000
#define T  4.0
#define A  2.0
#define M  121
#define TAU 1.0
#define DEC 50

int main(int argc, char* argv[])
{

	double*    t  = NULL;    // время (сек)
	double*    s  = NULL;    // входной сигнал
	complex_t* sc = NULL;
	complex_t S[M];    // комплексный спектр периодического сигнала
	double    Smag[M]; // амплитудный спектр периодического сигнала
	double	  Sphi[M]; // фазовый спектр периодического сигнала
	double    w[M];    // частота (рад/c) дискретного спектра

	void*	  handle;
	int       k, i;
	handle = dspl_load();
	if(!handle)
	{
		printf("cannot to load libdspl!\n");
		return 0;
	}

	t = (double*)malloc(N*sizeof(double));
	s = (double*)malloc(N*sizeof(double));
	sc = (complex_t*)malloc(N*sizeof(complex_t));

	// массив времени
	linspace(-T*2.0, T*2.0,  N, DSPL_PERIODIC, t);

        // 4 периода повторения п-импульса
        signal_pimp(t, N, A, TAU, 0.0, T, s);

	// расчет спектра
	fourier_series_dec(t, s, N, T, M, w, S);

	for(k = 0; k < M; k++)
	{
		Smag[k] = ABS(S[k])/4.0;
		Sphi[k] = atan2(IM(S[k]), RE(S[k]));
	}
	writetxt(w, Smag, M, "dat/fourier_series_prop_spectrum_amp_a.txt");
	writetxt(w, Sphi, M, "dat/fourier_series_prop_spectrum_phi_a.txt");

	double w0 = 10.0*M_PI;
	for(k = 0; k < N; k++)
	{
		RE(sc[k]) = s[k] * cos(w0*t[k]);
		IM(sc[k]) = s[k] * sin(w0*t[k]);
	}



	// разложение в ряд Фурье
        fourier_series_dec_cmplx(t, sc, N, T, M, w, S);
        for(k = 0; k < M; k++)
	{
		Smag[k] = ABS(S[k])/4.0;
		Sphi[k] = atan2(IM(S[k]), RE(S[k]));
	}
	writetxt(w, Smag, M, "dat/fourier_series_prop_spectrum_amp_se.txt");
	writetxt(w, Sphi, M, "dat/fourier_series_prop_spectrum_phi_se.txt");




	// сохраняю в файл dat/fourier_series_prop_time_a.txt
	decimate(t, N, DEC, t, &i);
	decimate(s, N, DEC, s, &i);
	writetxt(t, s, i, "dat/fourier_series_prop_time_a.txt");

	decimate_cmplx(sc, N, DEC, sc, &i);
	writetxt_cmplx_re(t, sc, i, "dat/fourier_series_prop_time_se_re.txt");
	writetxt_cmplx_im(t, sc, i, "dat/fourier_series_prop_time_se_im.txt");

	// массив времени
	linspace(-T*2.0, T*2.0,  N, DSPL_PERIODIC, t);
	// 4 периода повторения п-импульса
        signal_pimp(t, N, A, TAU, 0.0, T, s);

	for(k = 0; k < N; k++)
	{
		s[k] *= cos(w0*t[k]);
	}


	// разложение в ряд Фурье
        fourier_series_dec(t, s, N, T, M, w, S);
        for(k = 0; k < M; k++)
	{
		Smag[k] = ABS(S[k])/4.0;
		Sphi[k] = atan2(IM(S[k]), RE(S[k]));
	}
	writetxt(w, Smag, M, "dat/fourier_series_prop_spectrum_amp_sc.txt");
	writetxt(w, Sphi, M, "dat/fourier_series_prop_spectrum_phi_sc.txt");


	// сохраняю в файл dat/fourier_series_prop_time_sc.txt
	decimate(t, N, DEC, t, &i);
	decimate(s, N, DEC, s, &i);
	writetxt(t, s, i, "dat/fourier_series_prop_time_sc.txt");


	free(s);
	free(t);
	free(sc);
	return 0;
}
