#include <stdio.h>
#include <string.h>
#include "dspl.h"

#define N  30
#define FS 120
#define F0 -20

int main(int argc, char* argv[]) 
{

    void* handle;
	handle = dspl_load();

    complex_t s[N]; // Входной сигнал
    complex_t X[N]; // ДПФ
    double f[N];    // Индексы спектральных отсчетов
    double S[N];    // Амплитудный спектр без перестановки
    double Ssh[N];  // Амплитудный спектр после перестановки
    int n;

    // входной сигнал 
    for(n = 0; n < N; n++)
    {
        RE(s[n]) = cos(M_2PI * (double)n * F0 / FS);
        IM(s[n]) = sin(M_2PI * (double)n * F0 / FS);
    }

    // ДПФ
    dft_cmplx(s, N, X);
    
    // Амплитудный спектр
    for(n = 0; n < N; n++)
        S[n] = ABS(X[n]);
   
    // Перестановка спектральных отсчетов
    fft_shift(S, N, Ssh);

    // заполнение массива индексов спектральных отсчетов
    linspace(0, N,  N, DSPL_PERIODIC, f);

    //сохранить данные для построения графика
    writetxt(f, S,   N, "dat/dft_freq_fig1.txt");
    writetxt(f, Ssh, N, "dat/dft_freq_fig5.txt");

    // remember to free the resource
    dspl_free(handle);

    return 0;
}
