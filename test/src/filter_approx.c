#include <stdio.h>
#include <string.h>
#include <math.h>
#include "dspl.h"

#define N  8000



void r2hdb(double* r, double* h, double* hdb, double ep2, int n)
{
	int k;
	for(k = 0; k < n; k++)
	{
		r[k] *= r[k];
		h[k] = 1.0/(1.0 + ep2 * r[k]);
		hdb[k] = 10.0*log10(h[k]);
	}	
}


int main(int argc, char* argv[]) 
{
	double w[N];    // время (сек)    
	double r[N];    // входной сигнал
	double h[N], hdb[N];
	double ep2, Rp, Rs, es2, m;
	int ord, k;
	
	void* handle;           // DSPL handle
	handle = dspl_load();   // Load DSPL function
	
	// заполняю массив частот в логарифмическом формате от 0.01 до 100
	logspace(-2, 2,  N, DSPL_PERIODIC, w); 
	
	
	ord = 4;    // порядок фильтра
	Rp = 0.3;   // Неравномерность в полосе пропускания (дБ)
	Rs = 15.0;  // Уровень подавления в полосе заграждения (дБ)
	
	// Параметры ep^2 и es^2
	ep2 = pow(10.0, Rp*0.1) - 1.0;
	es2 = pow(10.0, Rs*0.1) - 1.0;
	
	// вывод на печать параметров ep^2 и es^2
	printf("ep^2 = %.4f\n", ep2);
	printf("es^2 = %.4f\n", es2);
	
	
	
	
	//**********************************************************************
	// Расчет F_N^2(w) и |H(jw)|^2 фильтра Баттерворта
	//**********************************************************************
	for(k = 0; k < N; k++)
	{
		r[k] = pow(w[k], (double)(ord));
	}
	r2hdb(r, h, hdb, ep2, N);
	
	// сохранение в файлы результатов расчета фильтра Баттерворта
	writetxt(w, r,   N, "dat/butter_r.txt");
	writetxt(w, h,   N, "dat/butter_h.txt");
	writetxt(w, hdb, N, "dat/butter_hdb.txt");
	


	
	
	//**********************************************************************
	// Расчет F_N^2(w) и |H(jw)|^2 фильтра Чебышева 1-го рода
	//**********************************************************************
	cheby_poly1(w, N, ord, r);
	r2hdb(r, h, hdb, ep2, N);
	// сохранение в файлы результатов расчета фильтра Чебышева 1-го рода    
	writetxt(w, r,   N, "dat/cheby1_r.txt");
	writetxt(w, h,   N, "dat/cheby1_h.txt");
	writetxt(w, hdb, N, "dat/cheby1_hdb.txt");
	
	
	
	
	
	//**********************************************************************
	 // Расчет F_N^2(w) и |H(jw)|^2 фильтра Чебышева 2-го рода
	//**********************************************************************
	for(k = 0; k < N; k++)
	{
		w[k] = 1.0 / w[k];
	}
	cheby_poly1(w, N, ord, r);
	for(k = 0; k < N; k++)
	{
		r[k]    =1.0 / (r[k] *r[k]);
		h[k]   = 1.0/(1.0 + es2*r[k]);
		hdb[k] = 10.0*log10(h[k]);
	}
	logspace(-2, 2,  N, DSPL_PERIODIC, w);
	// сохранение в файлы результатов расчета фильтра Чебышева 1-го рода
	writetxt(w, r,   N, "dat/cheby2_r.txt");
	writetxt(w, h,   N, "dat/cheby2_h.txt");
	writetxt(w, hdb, N, "dat/cheby2_hdb.txt");
	
	
	
	
	
	//**********************************************************************
	// Расчет F_N^2(w) и |H(jw)|^2 эллиптического фильтра
	//**********************************************************************
	ellip_modulareq(Rp, Rs, ord, &m);   // пересчет эллиптического модуля
	printf("modular m = %.3f\n", m);    // вывод на печать
	// расчет эллиптической рациональной функции
	ellip_rat(w, N, ord, m, r);
	r2hdb(r, h, hdb, ep2, N);
	// сохранение в файлы результатов расчета эллиптического фильтра   
	writetxt(w, r, N, "dat/ellip_r.txt");
	writetxt(w, h, N, "dat/ellip_h.txt");
	writetxt(w, hdb, N, "dat/ellip_hdb.txt");

	
	
	dspl_free(handle);      // free dspl handle
	return 0;
}



