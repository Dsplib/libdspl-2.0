/*
* Copyright (c) 2015-2019 Sergey Bakhurin
* Digital Signal Processing Library [http://dsplib.org]
*
* This file is part of DSPL.
*
* is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* DSPL is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with Foobar.    If not, see <http://www.gnu.org/licenses/>.
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"


#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int verif(double* x,    double* y, size_t n, double eps, double* err)
\brief Real arrays verification

Function calculates a maximum relative error between two real arrays `x`
and `y` (both length equals `n`):

\f[
e = \max \left( \frac{|x(k) - y(k)| }{ |x(k)|} \right),
\quad if \quad |x(k)| > 0,
\f]
or
\f[
e = \max(|x(k) - y(k)| ), ~\qquad if \quad~|x(k)| = 0,
\f]

This function can be used for algorithms verification if vector `x` is user
algorithm result and vector `y` -- reference vector.

\param[in] x
Pointer to the first vector `x`. \n
Vector size is `[n x 1]`. \n \n

\param[in] y
Pointer to the second vector `y`. \n
Vector size is `[n x 1]`. \n \n

\param[in] n
Size of vectors `x` and `y`. \n \n

\param[in] eps
Relative error threshold. \n
If error less than `eps`, then function returns
`DSPL_VERIF_SUCCESS`, else `DSPL_VERIF_FAILED`.    \n \n

\param[in, out] err
Pointer to the variable which keep
maximum relative error. \n
Pointer can be `NULL`, maximum error will not be returned
in this case. \n \n

\return
`DSPL_VERIF_SUCCESS` if maximum relative error less than `eps`. \n
Otherwise `DSPL_VERIF_FAILED`.

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int verif(double* x, double* y, size_t n, double eps, double* err)
\brief Верификация вещественных массивов

Функция производит расчет максимальной относительной ошибки между вещественными
векторами `x` и `y` равной длины `n`:

\f[
e = \max \left( \frac{|x(k) - y(k)| }{ |x(k)|} \right), \quad    \quad |x(k)| > 0,
\f]
или
\f[
e = \max(|x(k) - y(k)| ), \qquad \quad~|x(k)| = 0,
\f]
и возвращает `DSPL_VERIF_SUCCESS` если
разница \f$ e\f$ меньше `eps`.
В противном случае возвращает `DSPL_VERIF_FAILED`. \n
Данная функция используется для верификации работы алгоритмов если вектор `x`
результат работы алгоритма пользователя, а `y` -- результат работы этого же
алгоритма сторонней функцией.

\param[in] x
Указатель на первый вектор `x`. \n
Размер вектора `[n x 1]`. \n \n

\param[in] y
Указатель на второй вектор `y`. \n
Размер вектора `[n x 1]`. \n \n

\param[in] n
Размер векторов `x` и `y`. \n \n

\param[in] eps
Допустимая относительная ошибка. \n
Если максимальная относительная ошибка меньше `eps`, то функция возвращает
`DSPL_VERIF_SUCCESS`, в противном случае `DSPL_VERIF_FAILED`.    \n \n

\param[in, out] err
Указатель на переменную максимальной относительной ошибки. \n
По данному адресу будет записано значение максимальной относительной ошибки. \n
Указатель может быть `NULL`, значение ошибки в этом случае
не возвращается. \n \n

\return
`DSPL_VERIF_SUCCESS` если относительная ошибка меньше `eps`. \n
 В противном случае `DSPL_VERIF_FAILED`.

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API verif(double* x,    double* y, size_t n, double eps, double* err)
{
    double d, maxd;
    size_t k;
    int res;
    if(!x || !y)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;
    if(eps <= 0.0 )
        return ERROR_NEGATIVE;

    maxd = -100.0;

    for(k = 0; k < n; k++)
    {
        d = fabs(x[k] - y[k]);
        if(fabs(x[k]) > 0.0)
        {
            d = d / fabs(x[k]);
            if(d > maxd)
                maxd = d;
        }
    }
    if(err)
        *err = maxd;

    if(maxd > eps)
        res = DSPL_VERIF_FAILED;
    else
        res = DSPL_VERIF_SUCCESS;

    return res;
}





#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API verif_data_gen(int len, int type, char* fn)
{
    double    *pd = NULL;
    complex_t *pc = NULL;
    random_t rnd = {0};
    int err;
    if(len < 1)
        return ERROR_SIZE;
    if(!fn)
        return ERROR_FNAME;
    
    err = random_init(&rnd, RAND_TYPE_MRG32K3A, NULL);
    if(err != RES_OK)
        goto exit_label;
    switch(type & DAT_MASK)
    {
        case DAT_DOUBLE:
            pd = (double*)malloc(len*sizeof(double));
            if(!pd)
            {
                err = ERROR_MALLOC;
                goto exit_label;
            }
            err  = randn(pd, len, 1.0, 10.0, &rnd);
            if(err != RES_OK)
                goto exit_label;
            
            err = writebin(pd, len, type, fn);
            if(err != RES_OK)
                goto exit_label;
            break;
        case DAT_COMPLEX:
            pc = (complex_t*)malloc(len*sizeof(complex_t));
            if(!pc)
            {
                err = ERROR_MALLOC;
                goto exit_label;
            }
            err  = randn((double*)pc, 2*len, 1.0, 10.0, &rnd);
            if(err != RES_OK)
                goto exit_label;
            
            err = writebin(pc, len, type, fn);
            if(err != RES_OK)
                goto exit_label;
            break;
        default:
            err = ERROR_DAT_TYPE;
    }
    
exit_label:
    if(pd)
        free(pd);
    if(pc)
        free(pc);
    return err;
}





#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int verif_cmplx(complex_t* x, complex_t* y, size_t n,
                    double eps, double* err)
\brief
Complex arrays verification

Function calculates a maximum relative error between two complex arrays `x`
and `y` (both length equals `n`):

\f[
e = \max \left( \frac{|x(k) - y(k)| }{ |x(k)|} \right),
\quad if \quad |x(k)| > 0,
\f]
or
\f[
e = \max(|x(k) - y(k)| ), ~\qquad if \quad~|x(k)| = 0,
\f]
Return `DSPL_VERIF_SUCCESS` if maximum relative error \f$ e\f$ less than `eps`.
Else returns `DSPL_VERIF_FAILED`. \n

This function can be used for algorithms verification if vector `x` is user
algorithm result and vector `y` -- reference vector.

\param[in] x
Pointer to the first vector `x`. \n
Vector size is `[n x 1]`. \n \n

\param[in] y
Pointer to the second vector `y`. \n
Vector size is `[n x 1]`. \n \n

\param[in] n
Size of vectors `x` and `y`. \n \n

\param[in] eps
Relative error threshold. \n
If error less than `eps`, then function returns
`DSPL_VERIF_SUCCESS`, else `DSPL_VERIF_FAILED`.    \n \n

\param[in, out] err
Pointer to the variable which keep
maximum relative error. \n
Pointer can be `NULL`, maximum error will not be returned
in this case. \n \n

\return
`DSPL_VERIF_SUCCESS` if maximum relative error less than `eps`. \n
Otherwise `DSPL_VERIF_FAILED`.

\author
Sergey Bakhurin
www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int verif_cmplx(complex_t* x,    complex_t* y, size_t n,
                    double eps, double* err)
\brief Верификация комплексных массивов

Функция производит расчет максимальной относительной ошибки между комплексными
векторами `x` и `y` равной длины `n`:

\f[
e = \max \left( \frac{|x(k) - y(k)|}{|x(k)|} \right), \quad    \quad |x(k)| > 0,
\f]
или
\f[
e = \max(|x(k) - y(k)| ), ~\qquad    \quad~|x(k)| = 0,
\f]
и возвращает `DSPL_VERIF_SUCCESS` если
разница \f$ e\f$ меньше `eps`.
В противном случае возвращает `DSPL_VERIF_FAILED`. \n
Данная функция используется для верификации работы алгоритмов если вектор `x`
результат работы алгоритма пользователя, а `y` -- результат работы этого же
алгоритма сторонней функцией.

\param[in] x
Указатель на первый вектор `x`. \n
Размер вектора `[n x 1]`. \n \n

\param[in] y
Указатель на второй вектор `y`. \n
Размер вектора `[n x 1]`. \n \n

\param[in] n
Размер векторов `x` и `y`. \n \n

\param[in] eps
Допустимая относительная ошибка. \n
Если максимальная относительная ошибка меньше `eps`, то    функция возвращает
`DSPL_VERIF_SUCCESS`, в противном случае `DSPL_VERIF_FAILED`.    \n \n

\param[in, out] err
Указатель на переменную максимальной относительной ошибки. \n
По данному адресу будет записано значение максимальной относительной ошибки. \n
Указатель может быть `NULL`, значение ошибки в этом
случае не возвращается. \n \n

\return
`DSPL_VERIF_SUCCESS` если функция выполнена успешно. \n
 В противном случае `DSPL_VERIF_FAILED`.

\author
Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API verif_cmplx(complex_t* x,    complex_t* y, size_t n,
                 double eps, double* err)
{

    complex_t d;
    double mx, md, maxd;
    size_t k;
    int res;
    if(!x || !y)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;
    if(eps <= 0.0 )
        return ERROR_NEGATIVE;

    maxd = -100.0;

    for(k = 0; k < n; k++)
    {
        RE(d) = RE(x[k]) - RE(y[k]);
        IM(d) = IM(x[k]) - IM(y[k]);
        md = ABS(d);
        mx = ABS(x[k]);
        if(mx > 0.0)
        {
            md = md / mx;
            if(md > maxd)
                maxd = md;
        }
    }
    if(err)
        *err = maxd;

    if(maxd > eps)
        res = DSPL_VERIF_FAILED;
    else
        res = DSPL_VERIF_SUCCESS;

    return res;
}




#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
void DSPL_API verif_str(double* yout, int nout, 
                        char* str_msg, char* outfn, char* logfn)
{
    char str[VERIF_STR_BUF] = {0};
    char msg[VERIF_STR_BUF] = {0};
    double *y = NULL;
    double derr = 0.0;
    int n,  verr, type;

    sprintf(str, "%s", str_msg);
    while(strlen(str) < VERIF_STR_LEN)
        str[strlen(str)] = VERIF_CHAR_POINT;

    readbin(outfn, (void**)(&y), &n, &type);

    if(nout != n)
    {
        sprintf(msg, "FAILED (out size [%d] neq [%d])", n, nout);
        strcat(str, msg);
        addlog(str, logfn);
        printf("%s\n", str);
        return;
    }
    
    if(type!=DAT_DOUBLE)
    {
        sprintf(msg, "FAILED (type is complex)");
        strcat(str, msg);
        addlog(str, logfn);
        printf("%s\n", str);
        return;
    }

    verr = verif(yout, y, nout, VERIF_LEVEL_DOUBLE, &derr);
    if(verr == DSPL_VERIF_SUCCESS)
        sprintf(msg, "ok (err = %12.4E)", derr);
    else
        sprintf(msg, "FAILED (err = %12.4E)", derr);
    strcat(str, msg);
    addlog(str, logfn);
    printf("%s\n", str);
}







#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
void DSPL_API verif_str_cmplx(complex_t* yout, int nout, 
                              char* str_msg, char* outfn, char* logfn)
{
    char str[VERIF_STR_BUF] = {0};
    char msg[VERIF_STR_BUF] = {0};
    complex_t *y = NULL;
    double derr = 0.0;
    int n,  verr, type;

    sprintf(str, "%s", str_msg);
    while(strlen(str) < VERIF_STR_LEN)
        str[strlen(str)] = VERIF_CHAR_POINT;

    readbin(outfn, (void**)(&y), &n, &type);

    if(nout != n)
    {
        sprintf(msg, "FAILED (out size [%d] neq [%d])", n, nout);
        strcat(str, msg);
        addlog(str, logfn);
        printf("%s\n", str);
        return;
    }
    
    if(type!=DAT_COMPLEX)
    {
        sprintf(msg, "FAILED (type is complex)");
        strcat(str, msg);
        addlog(str, logfn);
        printf("%s\n", str);
        return;
    }

    verr = verif_cmplx(yout, y, nout, VERIF_LEVEL_DOUBLE, &derr);
    if(verr == DSPL_VERIF_SUCCESS)
        sprintf(msg, "ok (err = %12.4E)", derr);
    else
        sprintf(msg, "FAILED (err = %12.4E)", derr);
    strcat(str, msg);
    addlog(str, logfn);
    printf("%s\n", str);
}



