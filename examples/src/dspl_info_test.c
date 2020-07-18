#include <stdio.h>
#include <stdlib.h>
#include "dspl.h"


int main(int argc, char* argv[])
{
    /* libdspl handle                                      */
    void* hdspl;

    /* Load libdspl functions                              */
    hdspl = dspl_load();

    /* Check libdspl handle.                               */
    /* If hdspl == NULL means problem with libdspl loading */
    if(!hdspl)
    {
        printf("libdspl loading error!\n");
        return -1;
    }

    /* Print libdspl info                                  */
    dspl_info();


    /* free dspl handle                                    */
    dspl_free(hdspl);
    return 0;
}

