#include <stdio.h>
#include <string.h>
#include <math.h>
#include "dspl.h"

#define N  1000

int main(int argc, char* argv[]) 
{

    double w[N], r[N], h[N];      
    double ep2, Gp2;
    int ord, k;
    void* handle;
	handle = dspl_load();
    if(!handle)
    {
        printf("cannot to load libdspl!\n");
        return 0;
    }    


    linspace(0, 2.5,  N, DSPL_PERIODIC, w);

    ord = 3;
    Gp2 = 0.9;
    ep2 = 1.0 / Gp2 -1; 



    cheby_poly1(w, N, ord, r);
    for(k = 0; k < N; k++)
    {
        r[k] *= r[k];
        h[k] = 6.0/(1.0 + ep2*r[k]);
        w[k] *= 4; 
    }

        
    writetxt(w, h, N, "dat/cheby1_approx.txt");   
    

    // remember to free the resource
    dspl_free(handle);

    return 0;
}
