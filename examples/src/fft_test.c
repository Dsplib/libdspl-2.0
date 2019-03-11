#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N 14

int main()
{
  void* handle;           // DSPL handle
  handle = dspl_load();   // Загрузка DSPL 

  double    x[N];         // массив входного сигнала
  complex_t y[N];         // массив результата БПФ
  fft_t pfft;             // FFT объект
  
  memset(&pfft, 0, sizeof(fft_t)); // Заполняем FFT структуру нулями
  
  fft_create(&pfft, N);            // Создаем FFT структуру для длины N
  
  // заполняем массив входного сигнала
  for(int k = 0; k < N; k++)
    x[k] = (double)k;
  
  fft(x, N, &pfft, y);            // FFT

  // Печать результата
  for(int k = 0; k < N; k++)
    printf("y[%2d] = %9.3f%9.3f\n", k, RE(y[k]), IM(y[k]));

  fft_free(&pfft);        // Очищаем структуру fft_t
  dspl_free(handle);      // Очищаем dspl handle
  return 0;
}


