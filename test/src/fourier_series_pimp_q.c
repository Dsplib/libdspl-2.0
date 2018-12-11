#include <stdio.h>
#include <string.h>
//#include "common.h"
#include "dspl.h"
//#include "plot.h"

#define N  1000
#define T  4.0
#define A  2.0
#define M  41

int main(int argc, char* argv[])
{

    double    t[N];    // время (сек)
    double    s[N];    // входной сигнал
    complex_t S[M];    // комплексный спектр периодического сигнала
    double    Smag[M]; // амплитудный спектр периодического сигнала
    double    w[M];    // частота (рад/c) дискретного спектра
    double    wc[N];   // частота (рад/с) огибающей спектра
    double    Sc[N];   // огибающая спектра
    double    tau;     // длительность импульса



    // скважность
    double    Q[3] = {5.0, 2.0, 1.25};

    int q, m, n;
    char fname[64];
    
	void* handle;
	
  handle = dspl_load();
	if(!handle)
	{
		printf("cannot to load libdspl!\n");
		return 0;
	}


    for(q = 0; q < 3; q++)
    {
        tau = T/Q[q];

        // заполнение массива временных отсчетов
        // на 4-x периодах повторения сигнала
        // для отображения на осциллограмме
        linspace(-T*2.0, T*2.0,  N, DSPL_PERIODIC, t);

        // 4 периода повторения п-импульса скважности Q[q]
        signal_pimp(t, N, A, tau, 0.0, T, s);

        // сохранение в текстовый файл
        sprintf(fname, "dat/pimp_time_%.2lf.csv", Q[q]);
        writetxt(t, s, N, fname);


        // заполнение массива временных отсчетов
        // на одном периоде повторения сигнала
        linspace(-T/2.0, T/2.0,  N, DSPL_PERIODIC, t);

        // один период повторения п-импульса скважности Q[q]
        signal_pimp(t, N, A, tau, 0.0, T, s);

        // разложение в ряд Фурье
        fourier_series_dec(t, s, N, T, M, w, S);

        // Рассчет амплитудного спектра
        for(m = 0; m < M; m++)
        {
            printf("S[%d] = %f     %f\n", m, RE(S[m]), IM(S[m]));
            Smag[m] = ABS(S[m]);
        }
        // Сохранение в файл амплитудного спетра для скважности Q[q]
        sprintf(fname, "dat/pimp_freq_discrete_%.2lf.csv", Q[q]);
        writetxt(w, Smag, M, fname);

        // Вектор частот непрерывной огибаюхей вида sin(w/2*tau) / (w/2*T)
        linspace(w[0], w[M-1],  N, DSPL_SYMMETRIC, wc);

        // Расчет огибающей
        for(n = 0; n < N; n++)
            Sc[n] = (wc[n] == 0.0) ? A/Q[q] : fabs( A * sin(0.5*wc[n]*tau) / (0.5*wc[n] * T));

        // сохранение огибающей в файл для скважности Q[q]
        sprintf(fname, "dat/pimp_freq_cont_%.2lf.csv", Q[q]);
        writetxt(wc, Sc, N, fname);
    }
	// remember to free the resource
	dspl_free(handle);
    return 0;
}
