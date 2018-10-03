#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dspl.h"

// размер массивов исходных и восстановленных сигналов
#define N  1000
// Период
#define T  2
// Длительность импульса
#define TAU 1
// Максимальное количество коэффициентов ряда Фурье
#define MAX_M 61
// Количество периодов
#define P 3.0

int main(int argc, char* argv[])
{
	double t[N];	// время
	double s[N];	// исходный сигнал
	double x[N];	// восстановленный сигнал
	double w[MAX_M];// массив частоты
	complex_t S[MAX_M];	// Спектр

	complex_t xc[N];	// восстановленный по спектру
				// комплексный сигнал

	// количество спектральных гармоник усеченного ряда
	// Заметим, что 5 гармоник усеченного ряда содержат
	// две пары комплексно-сопряженных спектральных компонент
	// и одну постоянную составляющую.
	// Это означает, что ряд из 5 гармоник в комплексной форме соответствует
	// двум значениям a_n и b_n ряда в тригонометрической форме.
	// Аналогично M = 21 соответствует 10 значениям a_n и b_n ряда
	// в тригонометрической форме.
	int M[4] = {5, 9, 21, MAX_M};

	void* handle;
	int k, n;
	char fname[64];	// имя файла


	// Загрузка libdspl
	handle = dspl_load();
	if(!handle)
	{
		printf("cannot to load libdspl!\n");
		return 0;
	}

	// Заполняю вектор времени
	linspace(-T*0.5*P, T*0.5*P, N, DSPL_PERIODIC, t);

	// последовательность прямоугольных импульсов
	signal_pimp(t, N, 1.0, TAU, 0.0, T, s);

	// сохраняем в файл dat/fourier_series_pimp0.txt
	writetxt(t, s, N, "dat/fourier_series_pimp0.txt");

	// цикл по разному количеству гармоник усеченного ряда
	for(k = 0; k < 4; k++)
	{
		// расчет спектра для текущего усеченного ряда
		fourier_series_dec(t, s, N, T, M[k], w, S);

		// нормировка потому что P периодов в исходном сигнале
		for(n = 0; n < M[k]; n++)
		{
			RE(S[n]) /= P;
			IM(S[n]) /= P;
		}

		// восстанавливаю сигнал по усеченному ряду
		fourier_series_rec(w, S, M[k], t, N, xc);

		// Комплексный восстановленный сигнал имеет очень маленькую
		// мнимую часть, поэтому просто берем в вектор x реальную часть
		cmplx2re(xc, N, x, NULL);

		// сохраняю в файл для последующего построения графика
		sprintf(fname, "dat/fourier_series_pimp_rec_%d.txt", M[k]);
		writetxt(t, x, N, fname);
	}


	// Пилообразный сигнал
	signal_saw(t, N, 0.5, 0.0, T, s);
	for(n = 0; n < N; n++)
		s[n] += 0.5;

	// сохраняем в файл dat/fourier_series_saw0.txt
	writetxt(t, s, N, "dat/fourier_series_saw0.txt");

	// цикл по разному количеству гармоник усеченного ряда
	for(k = 0; k < 4; k++)
	{
		// расчет спектра для текущего усеченного ряда
		fourier_series_dec(t, s, N, T, M[k], w, S);

		// нормировка потому что P периодов в исходном сигнале
		for(n = 0; n < M[k]; n++)
		{
			RE(S[n]) /= P;
			IM(S[n]) /= P;
		}

		// восстанавливаю сигнал по усеченному ряду
		fourier_series_rec(w, S, M[k], t, N, xc);

		// Комплексный восстановленный сигнал имеет очень маленькую
		// мнимую часть, поэтому просто берем в вектор x реальную часть
		cmplx2re(xc, N, x, NULL);

		// сохраняю в файл для последующего построения графика
		sprintf(fname, "dat/fourier_series_saw_rec_%d.txt", M[k]);
		writetxt(t, x, N, fname);
	}

	// remember to free the resource
	dspl_free(handle);

	return system("gnuplot -p  gnuplot/fourier_series_rec.plt");
}
