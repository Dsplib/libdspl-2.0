#include <stdio.h>
#include <stdlib.h>
#include "dspl.h"


int main()
{
    void* handle;           // DSPL handle
    handle = dspl_load();   // Load DSPL function
    
    dspl_info();            // Print     DSPL information
    
    dspl_free(handle);      // free dspl handle
	return 0;
}


