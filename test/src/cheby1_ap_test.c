#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

// Порядок фильтра
#define ORD 4

// размер векторов частотной характериситки фильтра
#define N   1000


int main()
{
	void* handle;           // DSPL handle
	handle = dspl_load();   // Load DSPL function
	
	double a[ORD+1], b[ORD+1];	// коэффицинеты H(s)
	double Rp = 3.0;	// неравномерность в полосе пропускания 3дБ
	// Частота (w), АЧХ (mag), ФЧХ (phi) и ГВЗ (tau)
	double w[N], mag[N], phi[N], tau[N];
	int k;

	// рассчитываем нормированный ФНЧ Чебышева 1 рода 
	int res = cheby1_ap(Rp, ORD, b, a);
	if(res != RES_OK)
		printf("error code = 0x%8x\n", res);
	
	// печать коэффициентов фильтра
	for(k = 0; k < ORD+1; k++)
		printf("b[%2d] = %9.3f     a[%2d] = %9.3f\n", k, b[k], k, a[k]);
	
	// вектор частоты в логарифмической шакале от 0.01 до 100 рад/c
	logspace(-2.0, 2.0, N , DSPL_SYMMETRIC, w);
	
	//частотные характеристика фильтра
	freqs_resp(b, a, ORD, w, N, DSPL_FLAG_LOG|DSPL_FLAG_UNWRAP, 
							mag, phi, tau);
	
	// Сохранить характеристики для построения графиков
	writetxt(w, mag, N, "dat/cheby1_ap_test_mag.txt");
	writetxt(w, phi, N, "dat/cheby1_ap_test_phi.txt");
	writetxt(w, tau, N, "dat/cheby1_ap_test_tau.txt");
	
	dspl_free(handle);      // free dspl handle
	
	// выполнить скрипт GNUPLOT для построения графиков 
	// по рассчитанным данным
	return system("gnuplot -p gnuplot/cheby1_ap_test.plt");;
}


