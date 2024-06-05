/*
* Copyright (c) 2015-2022 Sergey Bakhurin
* Digital Signal Processing Library [http://dsplib.org]
*
* This file is part of libdspl-2.0.
*
* is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser  General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* DSPL is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public License
* along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef DSPL_H
#define DSPL_H

#include <math.h>

/* math const definition */
#ifndef M_PI
    #define M_PI        3.1415926535897932384626433832795
#endif

#ifndef M_2PI
    #define M_2PI       6.283185307179586476925286766559
#endif


#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup TYPES_GROUP
\typedef complex_t
\brief Complex data type.

libdspl-2.0 describes complex numbers data type as an array
of two `double` elements.
First element sets real part, second --- imaginary part.

For example:

\code{.cpp}
    complex_t z;
    z[0] =  1.0;
    z[1] = -2.0;
\endcode

Variable `z = 1-2j`, here `j` - imaginary unit.

For the convenience of working with complex numbers implemented
special macros: \ref RE, \ref IM, \ref ABSSQR
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup TYPES_GROUP
\typedef complex_t
\brief Описание комплексного типа данных.

Комплексный тип данных в библиотеке libdspl-2.0 определен как
массив из двух элементов типа `double`.
При этом первый элемент массива определяет реальную часть
комплексного числа, а второй - мнимую.

Например:

\code{.cpp}
    complex_t z;
    z[0] =  1.0;
    z[1] = -2.0;
\endcode

Переменная `z = 1-2j`, где `j` - мнимая единица.

Для удобства работы с комплексными числами реализованы
специальные макросы: \ref RE, \ref IM, \ref ABSSQR
***************************************************************************** */
#endif
typedef double complex_t[2];



/* Point 2D point2d_t[0] - x
            point2d_t[1] - y
*/
typedef double  point2d_t[2];

typedef struct
{
    point2d_t* points; /* line points array */
    int npoints;       /* number of points  */
}line2d_t;

typedef struct
{
    line2d_t* lines;  /* lines array     */
    int nlines;       /* number of lines */
    double level;     /* contour level   */
}contour2d_t;




#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup DFT_GROUP
\struct fft_t
\brief Fast Fourier Transform Object Data Structure

The structure stores pointers to twiddle factors and arrays of intermediate
data of the fast Fourier transform algorithm.

The libdspl-2.0 library uses an FFT algorithm for composite size.

\param  n
The size of the FFT vector for which memory is allocated
in the structure arrays.  \n
The parameter `n` must be equal to an integer power of two (radix 2). \n \n

\param  w
Pointer to the vector of twiddle factors. \n
The size of the vector is `[n x 1]`. \n
The memory must be allocated and an array of twiddle factors
must be filled with the \ref fft_create function. \n\n

\param  t0
Pointer to the vector of intermediate results of the FFT algorithm. \n
The size of the vector is `[n x 1]`. \n
Memory must be allocated by \ref fft_create function. \n\n

\param  t1
Pointer to the vector of intermediate results. \n
The size of the vector is `[n x 1]`. \n
The memory must be allocated with the \ref fft_create function. \n\n

\param w32
Static twiddle factors vector for 32-points FFT. \n \n

\param w64
Static twiddle factors vector for 32-points FFT. \n \n

\param w128
Static twiddle factors vector for 32-points FFT. \n \n

\param w256
Static twiddle factors vector for 32-points FFT. \n \n

\param w512
Static twiddle factors vector for 32-points FFT. \n \n

\param w1024
Dynamic twiddle factors vector for 32-points FFT. \n \n

\param w2048
Dynamic twiddle factors vector for 32-points FFT. \n \n

\param w4096
Dynamic twiddle factors vector for 32-points FFT. \n \n

The structure is calculated with the \ref fft_create function once
before using the FFT algorithm. \n
A pointer to an object of this structure may be
reused when calling FFT functions. \n
Before exiting the program, dedicated memory for twiddle factors and arrays of
intermediate data must be cleared by the \ref fft_free function.

For example:

\code
fft_t pfft = {0};     // Structure fft_t and clear all fields
int n = 64;           // FFT size

int err;

// Create and fill FFT structure for 64-points FFT
err = fft_create(&pfft, n);

// FFT calculation here
// FFT calculation here one more
// ...

// Clear fft structure
fft_free(&pfft);
\endcode

\note
It is important to note that if the object `fft_t` was created for the FFT size
equal to` n`, it can only be used for FFT of size `n`. \n \n
It’s also worth noting that the FFT functions independently control the size,
and independently allocate the memory of the FFT object, if necessary.
So if you call any function using the `fft_t` structure with filled
data for the FFT length `k` for calculating the FFT of length`n`,
then the structure arrays will be automatically recreated for the length `n`.

\author  Sergey Bakhurin  www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup DFT_GROUP
\struct fft_t
\brief Структура данных объекта быстрого преобразования Фурье

Структура хранит указатели на массивы поворотных коэффициентов
и массивы промежуточных данных алгоритма быстрого преобразования Фурье.

Библиотека libdspl-2.0 использует для БПФ алгоритм для составной длины

\param  n
Размер вектора БПФ, для которого выделена память в массивах структуры.  \n
Парметр `n` должен быть равен целой степени двойки. \n \n

\param  w
Указатель на вектор поворотных коэффициентов алгоритма БПФ. \n
Размер вектора `[n x 1]`.  \n
Память должна быть выделена и массив поворотных коэффициентов
должен быть заполнен функцией \ref fft_create.  \n \n

\param  t0
Указатель на вектор промежуточных вычислений алгоритма БПФ. \n
Размер вектора `[n x 1]`. \n
Память должна быть выделена функцией \ref fft_create. \n \n

\param  t1
Указатель на вектор промежуточных вычислений алгоритма БПФ. \n
Размер вектора `[n x 1]`. \n
Память должна быть выделена функцией \ref fft_create. \n \n

\param w32
Статический вектор поворотных коэффициентов 32-точечного БПФ. \n \n

\param w64
Статический вектор поворотных коэффициентов 64-точечного БПФ. \n \n

\param w128
Статический вектор поворотных коэффициентов 128-точечного БПФ. \n \n

\param w256
Статический вектор поворотных коэффициентов 256-точечного БПФ. \n \n

\param w512
Статический вектор поворотных коэффициентов 512-точечного БПФ. \n \n

\param w1024
Статический вектор поворотных коэффициентов 1024-точечного БПФ. \n \n

\param w2048
Статический вектор поворотных коэффициентов 2048-точечного БПФ. \n \n

\param w4096
Статический вектор поворотных коэффициентов 4096-точечного БПФ. \n \n


Структура заполняется функцией \ref fft_create один раз
до использования алгоритма БПФ.  \n
Указатель на объект данной структуры может быть
многократно использован при вызове функций БПФ. \n
Перед выходом из программы выделенную память под поворотные
коэффициенты и массивы промежуточных данных
необходимо очистить функцией \ref fft_free. Например:
\code
fft_t pfft = {0};     // объявляем объект fft_t и обнуляем все поля
int n = 64;           // Размер БПФ
int err;

// создаем объект для 64-точечного БПФ
err = fft_create(&pfft, n);

// Вызов БПФ функции
// Еще раз вызов БПФ функции
// ...

// очистить память объекта БПФ
fft_free(&pfft);
\endcode

\note
Важно отметить, что если объект `fft_t` был создан для размера БПФ равного `n`,
то он может быть использован только для БПФ размера `n`.  \n\n
Также необходимо заметить, что функции БПФ самостоятельно контролируют размер,
и самостоятельно выделяют память объекта БПФ при необходимости.
Так если вызвать любую функцию использующую структуру `fft_t` с заполненными
данными для длины БПФ `k` для расчета БПФ длины `n`,
то массивы структуры будут автоматически пересозданы для длины `n`.

\author
Бахурин Сергей.
www.dsplib.org
***************************************************************************** */
#endif
typedef struct
{
    complex_t*  w;
    complex_t*  t0;
    complex_t*  t1;
   
    /* radix-2 twiddle factors vectors */
    complex_t    w32[ 32];
    complex_t    w64[ 64];
    complex_t   w128[128];
    complex_t   w256[256];
    complex_t   w512[512];
    complex_t*  w1024;
    complex_t*  w2048;
    complex_t*  w4096;
    int         n;
} fft_t;



#define RAND_TYPE_MRG32K3A 0x00000001
#define RAND_TYPE_MT19937  0x00000002
#define RAND_MT19937_NN    312

#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup SPEC_MATH_RAND_GEN_GROUP
\struct random_t
\brief Структура параметров датчиков псевдослучайных чисел.

Структура хранит инициализацию и текущие регистры различных датчиков
псевдослучайных чисел. В библиотеке используются следующие датчики:
\li MRG32K3A -- 32 битный датчик разработан Пьером Лекуэром [1].
\li MT19937-64 -- 64-битный датчик
<a href = "https://en.wikipedia.org/wiki/Mersenne_Twister">
Вихрь Мерсенна
</a> [2, 3].

\note
[1] Pierre L'Ecuyer, (1999) Good Parameters and Implementations for Combined
    Multiple Recursive Random Number Generators. Operations Research
    47(1):159-164. https://doi.org/10.1287/opre.47.1.159 \n\n
[2] T. Nishimura, ``Tables of 64-bit Mersenne Twisters // ACM Transactions
    on Modeling and Computer Simulation 10. (2000) 348--357. \n\n
[3] M. Matsumoto and T. Nishimura  Mersenne Twister: a 623-dimensionally
    equidistributed uniform pseudorandom number generator // ACM Transactions
    on Modeling and Computer Simulation 8. (Jan. 1998) 3--30.  \n\n

\param  mrg32k3a_seed
Начальная инициализация датчика MRG32K3A. \n \n

\param  mrg32k3a_x
Первый вектор состояния рекурсивного датчика MRG32K3A. \n \n

\param  mrg32k3a_y
Второй вектор состояния рекурсивного датчика MRG32K3A. \n \n

\param  mt19937_mt
Первый вектор состояния рекурсивного датчика MT19937-64. \n \n

\param  mt19937_mti
Текущий индекс в векторе состояния датчика MT19937-64. \n \n

Параметры данной структуры заполняются автоматически функцией `random_init`
и используются функциями генерации псевдослучайных векторов.

\author Бахурин Сергей. www.dsplib.org
***************************************************************************** */
#endif
typedef struct
{

    double mrg32k3a_seed;
    double mrg32k3a_x[3];
    double mrg32k3a_y[3];

    /* The array for the MT19937 state vector */
    unsigned long long mt19937_mt[RAND_MT19937_NN];
    int                mt19937_mti;

    int type;

}random_t;



#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup TYPES_GROUP
\def RE(x)
\brief Macro sets real part of the complex number.

Example:
\code{.cpp}
    complex_t z;
    RE(z) =  1.0;
    IM(z) = -2.0;
\endcode

Variable `z = 1-2j`, here `j` - imaginary unit.

This macro can be used to return
real part of the complex number:

\code{.cpp}
    complex_t z = {3.0, -4.0};
    double    r;
    r = RE(z);
\endcode
In this example `z = 3-4i`,
but variable `r` will keep 3.
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup TYPES_GROUP
\def RE(x)
\brief Макрос определяющий реальную часть комплексного числа.

Например:
\code{.cpp}
    complex_t z;
    RE(z) =  1.0;
    IM(z) = -2.0;
\endcode

Переменная `z = 1-2j`, где `j` --- мнимая единица.

Аналогично, макрос можно использовать для получения
реальной части комплексного числа:

\code{.cpp}
    complex_t z = {3.0, -4.0};
    double    r;
    r = RE(z);
\endcode
В данном примере переменная `z = 3-4i`, а в переменой `r`
будет храниться число 3.
***************************************************************************** */
#endif
#define RE(x) (x[0])



#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup TYPES_GROUP
\def IM(x)
\brief Macro sets imaginary part of the complex number.

Example:
\code{.cpp}
    complex_t z;
    RE(z) =  1.0;
    IM(z) = -2.0;
\endcode

Variable `z = 1-2j`, here `j` - imaginary unit.

This macro can be used to return
imaginary part of the complex number:
\code{.cpp}
    complex_t z = {3.0, -4.0};
    double    r;
    r = IM(z);
\endcode
In this example `z = 3-4i`,
but variable `r` will keep -4.
***************************************************************************** */


#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup TYPES_GROUP
\def IM(x)
\brief Макрос определяющий мнимую часть комплексного числа.

Например:
\code{.cpp}
    complex_t z;
    RE(z) =  1.0;
    IM(z) = -2.0;
\endcode

Переменная `z = 1-2j`, где `j` - мнимая единица.

Аналогично, макрос можно использовать для получения
мнимой части комплексного числа:
\code{.cpp}
    complex_t z = {3.0, -4.0};
    double r;
    r = IM(z);
\endcode
В данном примере переменная `z = 3-4i`,
а в переменой `r` будет храниться число -4.
***************************************************************************** */
#endif
#define IM(x) (x[1])



#define SQR(x) ((x) * (x))


#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup TYPES_GROUP
\def ABSSQR(x)
\brief
The macro returns the square of the modulus of a complex number `x`.

Square of the modulus of a complex number \f$ x = a + j  b \f$ equals:

\f[
    |x|^2 = x x^* = a^2 + b^2.
\f]

Example:
\code{.cpp}
    complex_t z;
    double y;
    RE(z) =  1.0;
    IM(z) = -2.0;
    y = ABSSQR(z);
\endcode

Variable `z = 1-2j`, here `j` - imaginary unit, but variable `y = 5`.
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup TYPES_GROUP
\def ABSSQR(x)
\brief Макрос возвращает квадрат модуля комплексного числа `x`.

Квадрат модуля комплексного числа \f$ x = a + j  b \f$ равен:

\f[
    |x|^2 = x x^* = a^2 + b^2.
\f]

Например:
\code{.cpp}
    complex_t z;
    double y;
    RE(z) =  1.0;
    IM(z) = -2.0;
    y = ABSSQR(z);
\endcode

Переменная `z = 1-2j`, где `j` - мнимая единица, а переменная `y = 5`.
***************************************************************************** */
#endif
#define ABSSQR(x) ((SQR(RE(x))) + (SQR(IM(x))))




#define ABS(x)          sqrt((ABSSQR(x)))


#define ARG(x)          atan2(IM(x), RE(x))


#define CMRE(a,b)       ((RE(a)) * (RE(b)) - (IM(a)) * (IM(b)))


#define CMIM(a,b)       ((RE(a)) * (IM(b)) + (IM(a)) * (RE(b)))


#define CMCONJRE(a, b)  ((RE(a)) * (RE(b)) + (IM(a)) * (IM(b)))


#define CMCONJIM(a, b)  ((IM(a)) * (RE(b)) - (RE(a)) * (IM(b)))



#define RES_OK                                0

/* Error codes                                          */
/* A                                          0x01xxxxxx*/
#define ERROR_ARG_PARAM                       0x01180716
/* B                                          0x02xxxxxx*/
/* C                                          0x03xxxxxx*/
/* D                                          0x04xxxxxx*/
#define ERROR_DAT_TYPE                        0x04012020
#define ERROR_DIV_ZERO                        0x04102226
/* E                                          0x05xxxxxx*/
#define ERROR_ELLIP_MODULE                    0x05121315
/* F                                          0x06xxxxxx*/
#define ERROR_FFT_SIZE                        0x06062021
#define ERROR_FILTER_A0                       0x06090100
#define ERROR_FILTER_APPROX                   0x06090116
#define ERROR_FILTER_FT                       0x06090620
#define ERROR_FILTER_ORD                      0x06091518
#define ERROR_FILTER_ORD_BP                   0x06091519
#define ERROR_FILTER_RP                       0x06091816
#define ERROR_FILTER_RS                       0x06091819
#define ERROR_FILTER_TYPE                     0x06092025
#define ERROR_FILTER_WP                       0x06092316
#define ERROR_FILTER_WS                       0x06092319
#define ERROR_FNAME                           0x06140113
#define ERROR_FOPEN                           0x06151605
#define ERROR_FREAD_SIZE                      0x06180501
#define ERROR_FS                              0x06190000
#define ERROR_FWRITE_SIZE                     0x06231820
/* G                                          0x07xxxxxx*/
#define ERROR_GNUPLOT_CREATE                  0x07161203
#define ERROR_GNUPLOT_FNPNG                   0x07161206
#define ERROR_GNUPLOT_TERM                    0x07161220
/* H                                          0x08xxxxxx*/
/* I                                          0x09xxxxxx*/
#define ERROR_INF                             0x09140600
/* J                                          0x10xxxxxx*/
/* K                                          0x11xxxxxx*/
/* L                                          0x12xxxxxx*/
#define ERROR_LAPACK                          0x12011601
/* M                                          0x13xxxxxx*/
#define ERROR_MALLOC                          0x13011212
#define ERROR_MATRIX_SIZE                     0x13011926
#define ERROR_MIN_MAX                         0x13091413
/* N                                          0x14xxxxxx*/
#define ERROR_NAN                             0x14011400
#define ERROR_NEGATIVE                        0x14050701
/* O                                          0x15xxxxxx*/
#define ERROR_OVERLAP                         0x15220412
/* P                                          0x16xxxxxx*/
#define ERROR_POLY_AN                         0x16150114
#define ERROR_POLY_ORD                        0x16151518
#define ERROR_PTR                             0x16201800
/* Q                                          0x17xxxxxx*/
/* R                                          0x18xxxxxx*/
#define ERROR_RAND_SIGMA                      0x18011909
#define ERROR_RAND_TYPE                       0x18012009
#define ERROR_RESAMPLE_RATIO                  0x18051801
#define ERROR_RESAMPLE_FRAC_DELAY             0x18050604
/* S                                          0x19xxxxxx*/
#define ERROR_SIZE                            0x19092605
#define ERROR_SYM_TYPE                        0x19251320
/* T                                          0x20xxxxxx*/
/* U                                          0x21xxxxxx*/
#define ERROR_UNWRAP                          0x21142318
/* V                                          0x22xxxxxx*/
/* W                                          0x23xxxxxx*/
#define ERROR_WIN_PARAM                       0x23091601
#define ERROR_WIN_SYM                         0x23091925
#define ERROR_WIN_TYPE                        0x23092025
/* X                                          0x24xxxxxx*/
#define ERROR_XCORR_FLAG                      0x24031518
/* Y                                          0x25xxxxxx*/
/* Z                                          0x26xxxxxx*/

#define DAT_MASK                              0x00000001
#define DAT_DOUBLE                            0x00000000
#define DAT_COMPLEX                           0x00000001

#define DSPL_MATRIX_BLOCK                     32


#define DSPL_SYMMETRIC                        0x00000000
#define DSPL_PERIODIC                         0x00000001

#define DSPL_FLAG_DIGITAL                     0x00000000
#define DSPL_FLAG_ANALOG                      0x00000001
#define DSPL_FLAG_LOGMAG                      0x00000002
#define DSPL_FLAG_UNWRAP                      0x00000004
#define DSPL_FLAG_FFT_SHIFT                   0x00000008
#define DSPL_FLAG_PSD_TWOSIDED                DSPL_FLAG_FFT_SHIFT




#define DSPL_WIN_SYM_MASK                     0x00000001
#define DSPL_WIN_MASK                         0x00FFFFFE

#define DSPL_WIN_SYMMETRIC                    DSPL_SYMMETRIC
#define DSPL_WIN_PERIODIC                     DSPL_PERIODIC


#define DSPL_WIN_BARTLETT                     0x00000004
#define DSPL_WIN_BARTLETT_HANN                0x00000008
#define DSPL_WIN_BLACKMAN                     0x00000010
#define DSPL_WIN_BLACKMAN_HARRIS              0x00000040
#define DSPL_WIN_BLACKMAN_NUTTALL             0x00000080
#define DSPL_WIN_FLAT_TOP                     0x00000100
#define DSPL_WIN_GAUSSIAN                     0x00000400
#define DSPL_WIN_HAMMING                      0x00000800
#define DSPL_WIN_HANN                         0x00001000
#define DSPL_WIN_LANCZOS                      0x00004000
#define DSPL_WIN_NUTTALL                      0x00008000
#define DSPL_WIN_RECT                         0x00010000
#define DSPL_WIN_COS                          0x00040000
#define DSPL_WIN_CHEBY                        0x00080000
#define DSPL_WIN_KAISER                       0x00100000


#define DSPL_FILTER_TYPE_MASK                 0x000000FF
#define DSPL_FILTER_LPF                       0x00000001
#define DSPL_FILTER_HPF                       0x00000002
#define DSPL_FILTER_BPASS                     0x00000004
#define DSPL_FILTER_BSTOP                     0x00000008

#define DSPL_FILTER_APPROX_MASK               0x0000FF00
#define DSPL_FILTER_BUTTER                    0x00000100
#define DSPL_FILTER_CHEBY1                    0x00000200
#define DSPL_FILTER_CHEBY2                    0x00000400
#define DSPL_FILTER_ELLIP                     0x00000800


#define DSPL_XCORR_NOSCALE                    0x00000000
#define DSPL_XCORR_BIASED                     0x00000001
#define DSPL_XCORR_UNBIASED                   0x00000002



#define ELLIP_ITER                            16
#define ELLIP_MAX_ORD                         24

#define  DSPL_VERIF_FAILED                    1
#define  DSPL_VERIF_SUCCESS                   0

#define PLOT_HOLD                             0x00000001

#define VERIF_STR_BUF         128
#define VERIF_STR_LEN         48
#define VERIF_CHAR_POINT      46
#define VERIF_LEVEL_COMPLEX   1E-11
#define VERIF_LEVEL_DOUBLE    1E-12



#ifdef __cplusplus
  extern "C" {
#endif



#ifdef BUILD_LIB
  /* Declare DSPL_API for Windows OS */
  #ifdef WIN_OS
    #define DSPL_API __declspec(dllexport)
  #endif /* WIN_OS */
  /* Declare DSPL_API for LINUX OS */
  #ifdef LINUX_OS
    #define DSPL_API
  #endif /* LINUX_OS */
#endif /* BUILD_DLL */

#define COMMA ,


#ifdef BUILD_LIB
    #define DECLARE_FUNC(type, fn, param)\
                         type DSPL_API fn(param);
#endif

#ifndef BUILD_LIB
    #define DECLARE_FUNC( type, fn, param)\
                          typedef type (*p_##fn)(param);\
                          extern p_##fn   fn;

#endif



/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        acos_cmplx,                  complex_t*
                                                COMMA int
                                                COMMA complex_t*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        addlog,                      char*         str
                                                COMMA char*         fn);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        array_scale_lin,             double*       x
                                                COMMA int           n
                                                COMMA double        xmin
                                                COMMA double        xmax
                                                COMMA double        dx
                                                COMMA double        h
                                                COMMA double*       y);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        asin_cmplx,                  complex_t*
                                                COMMA int
                                                COMMA complex_t*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        bessel_i0,                   double*       x
                                                COMMA int           n
                                                COMMA double*       y);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        bilinear,                    double*       bs
                                                COMMA double*       as
                                                COMMA int           ord
                                                COMMA double*       bz
                                                COMMA double*       az);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        butter_ap,                   double
                                                COMMA int
                                                COMMA double*
                                                COMMA double*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        butter_ap_zp,                int
                                                COMMA double
                                                COMMA complex_t*
                                                COMMA int*
                                                COMMA complex_t*
                                                COMMA int*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        cheby_poly1,                 double*
                                                COMMA int
                                                COMMA int
                                                COMMA double*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        cheby_poly2,                 double*
                                                COMMA int
                                                COMMA int
                                                COMMA double*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        cheby1_ap,                   double
                                                COMMA int
                                                COMMA double*
                                                COMMA double*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        cheby1_ap_zp,                int
                                                COMMA double
                                                COMMA complex_t*
                                                COMMA int*
                                                COMMA complex_t*
                                                COMMA int*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        cheby2_ap,                   double           rs
                                                COMMA int              ord
                                                COMMA double*          b
                                                COMMA double*          a);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        cheby2_ap_wp1,               double           rp
                                                COMMA double           rs
                                                COMMA int              ord
                                                COMMA double*          b
                                                COMMA double*          a);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        cheby2_ap_zp,                int
                                                COMMA double
                                                COMMA complex_t*
                                                COMMA int*
                                                COMMA complex_t*
                                                COMMA int*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        cmplx2re,                    complex_t*
                                                COMMA int
                                                COMMA double*
                                                COMMA double*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        concat,                      void*
                                                COMMA size_t
                                                COMMA void*
                                                COMMA size_t
                                                COMMA void*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        contour2d,                   double* z
                                                COMMA double* x 
                                                COMMA double* y
                                                COMMA int     n
                                                COMMA int     m
                                                COMMA double  lev 
                                                COMMA contour2d_t* c);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(void,       contour2d_free,              contour2d_t* c);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        conv,                        double*
                                                COMMA int
                                                COMMA double*
                                                COMMA int
                                                COMMA double*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        conv_cmplx,                  complex_t*
                                                COMMA int
                                                COMMA complex_t*
                                                COMMA int
                                                COMMA complex_t*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        conv_fft,                    double*           a
                                                COMMA int               na
                                                COMMA double*           b
                                                COMMA int               nb
                                                COMMA fft_t*            pfft
                                                COMMA int               nfft
                                                COMMA double*           c);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        conv_fft_cmplx,              complex_t*        a
                                                COMMA int               na
                                                COMMA complex_t*        b
                                                COMMA int               nb
                                                COMMA fft_t*            pfft
                                                COMMA int               nfft
                                                COMMA complex_t*        c);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        cos_cmplx,                   complex_t*
                                                COMMA int
                                                COMMA complex_t*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        decimate,                    double*           x
                                                COMMA int               n
                                                COMMA int               d
                                                COMMA double*           y
                                                COMMA int*              cnt);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        decimate_cmplx,              complex_t*        x
                                                COMMA int               n
                                                COMMA int               d
                                                COMMA complex_t*        y
                                                COMMA int*              cnt);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        dft,                         double*
                                                COMMA int
                                                COMMA complex_t*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        dft_cmplx,                   complex_t*
                                                COMMA int
                                                COMMA complex_t*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(double,     dmod,                        double
                                                COMMA double);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(void,       dspl_info,                   void);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        ellip_acd,                   double*           w
                                                COMMA int               n
                                                COMMA double            k
                                                COMMA double*           u);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        ellip_acd_cmplx,             complex_t*        w
                                                COMMA int               n
                                                COMMA double            k
                                                COMMA complex_t*        u);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        ellip_ap,                    double            rp
                                                COMMA double            rs
                                                COMMA int               ord
                                                COMMA double*           b
                                                COMMA double*           a);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        ellip_ap_zp,                 int               ord
                                                COMMA double            rp
                                                COMMA double            rs
                                                COMMA complex_t*        z
                                                COMMA int*              nz
                                                COMMA complex_t*        p
                                                COMMA int*              np);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        ellip_asn,                   double*           w
                                                COMMA int               n
                                                COMMA double            k
                                                COMMA double*           u);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        ellip_asn_cmplx,             complex_t*        w
                                                COMMA int               n
                                                COMMA double            k
                                                COMMA complex_t*        u);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        ellip_cd,                    double*           u
                                                COMMA int               n
                                                COMMA double            k
                                                COMMA double*           y);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        ellip_cd_cmplx,              complex_t*        u
                                                COMMA int               n
                                                COMMA double            k
                                                COMMA complex_t*        y);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        ellip_landen,                double            k
                                                COMMA int               n
                                                COMMA double*           y);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        ellip_modulareq,             double            rp
                                                COMMA double            rs
                                                COMMA int               ord
                                                COMMA double*           k);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        ellip_rat,                   double*           w
                                                COMMA int               n
                                                COMMA int               ord
                                                COMMA double            k
                                                COMMA double*           u);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        ellip_sn,                    double*           u
                                                COMMA int               n
                                                COMMA double            k
                                                COMMA double*           y);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        ellip_sn_cmplx,              complex_t*        u
                                                COMMA int               n
                                                COMMA double            k
                                                COMMA complex_t*        y);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        farrow_lagrange,             double*
                                                COMMA int
                                                COMMA double
                                                COMMA double
                                                COMMA double
                                                COMMA double**
                                                COMMA int*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        farrow_spline,               double*
                                                COMMA int
                                                COMMA double
                                                COMMA double
                                                COMMA double
                                                COMMA double**
                                                COMMA int*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        fft,                         double*
                                                COMMA int
                                                COMMA fft_t*
                                                COMMA complex_t* );
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        fft_abs,                     double*          x
                                                COMMA int              n
                                                COMMA fft_t*           pfft
                                                COMMA double           fs
                                                COMMA int              flag
                                                COMMA double*          mag
                                                COMMA double*          freq);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        fft_abs_cmplx,               complex_t*       x
                                                COMMA int              n
                                                COMMA fft_t*           pfft
                                                COMMA double           fs
                                                COMMA int              flag
                                                COMMA double*          mag
                                                COMMA double*          freq);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        fft_cmplx,                   complex_t*
                                                COMMA int
                                                COMMA fft_t*
                                                COMMA complex_t* );
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        fft_create,                  fft_t*
                                                COMMA int);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(void,       fft_free,                    fft_t*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        fft_mag,                     double*          x
                                                COMMA int              n
                                                COMMA fft_t*           pfft
                                                COMMA double           fs
                                                COMMA int              flag
                                                COMMA double*          mag
                                                COMMA double*          freq);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        fft_mag_cmplx,               complex_t*       x
                                                COMMA int              n
                                                COMMA fft_t*           pfft
                                                COMMA double           fs
                                                COMMA int              flag
                                                COMMA double*          mag
                                                COMMA double*          freq);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        fft_shift,                   double*
                                                COMMA int n
                                                COMMA double*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        fft_shift_cmplx,             complex_t*
                                                COMMA int
                                                COMMA complex_t*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        filter_freq_resp,            double*          b
                                                COMMA double*          a
                                                COMMA int              ord
                                                COMMA double*          w
                                                COMMA int              n
                                                COMMA int              flag
                                                COMMA double*          mag
                                                COMMA double*          phi
                                                COMMA double*          tau);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        filter_iir,                  double*
                                                COMMA double*
                                                COMMA int
                                                COMMA double*
                                                COMMA int
                                                COMMA double*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(double,     filter_ws1,                  int              ord
                                                COMMA double           rp
                                                COMMA double           rs
                                                COMMA int              type);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        filter_zp2ab,                complex_t*
                                                COMMA int
                                                COMMA complex_t*
                                                COMMA int
                                                COMMA int
                                                COMMA double*
                                                COMMA double*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        find_max_abs,                double*        a
                                                COMMA int            n
                                                COMMA double*        m
                                                COMMA int*          ind);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        find_nearest,                double*        x
                                                COMMA int            n
                                                COMMA double         val
                                                COMMA int*           idx
                                                COMMA double*        dist);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        fir_linphase,                int            ord
                                                COMMA double         w0
                                                COMMA double         w1
                                                COMMA int            filter_type
                                                COMMA int            wintype
                                                COMMA double         winparam
                                                COMMA double*        h);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        flipip,                      double*
                                                COMMA int);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        flipip_cmplx,                complex_t*
                                                COMMA int);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        fourier_integral_cmplx,      double*         t
                                                COMMA complex_t*      s
                                                COMMA int             nt
                                                COMMA int             nw
                                                COMMA double*         w
                                                COMMA complex_t*      y);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        fourier_series_dec,          double*
                                                COMMA double*
                                                COMMA int
                                                COMMA double
                                                COMMA int
                                                COMMA double*
                                                COMMA complex_t*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        fourier_series_dec_cmplx,    double*         t
                                                COMMA complex_t*      s
                                                COMMA int             nt
                                                COMMA double          period
                                                COMMA int             nw
                                                COMMA double*         w
                                                COMMA complex_t*      y);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        fourier_series_rec,          double*
                                                COMMA complex_t*
                                                COMMA int
                                                COMMA double*
                                                COMMA int
                                                COMMA complex_t*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        freqs,                       double*
                                                COMMA double*
                                                COMMA int
                                                COMMA double*
                                                COMMA int
                                                COMMA complex_t*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        freqs_cmplx,                 double*          b
                                                COMMA double*          a
                                                COMMA int              ord
                                                COMMA complex_t*       s
                                                COMMA int              n
                                                COMMA complex_t*       h);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        freqs2time,                  double*
                                                COMMA double*
                                                COMMA int
                                                COMMA double
                                                COMMA int
                                                COMMA fft_t*
                                                COMMA double*
                                                COMMA double*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        freqz,                       double*
                                                COMMA double*
                                                COMMA int
                                                COMMA double*
                                                COMMA int
                                                COMMA complex_t*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(void,       gnuplot_close,               void*             h);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(void,       gnuplot_cmd,                 void*             h
                                                COMMA char*             cmd);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        gnuplot_create,              int               argc
                                                COMMA char*             argv[]
                                                COMMA int               w
                                                COMMA int               h
                                                COMMA char*             fn_png
                                                COMMA void**            hplot);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        gnuplot_open,                void**            hplot);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        goertzel,                    double*
                                                COMMA int
                                                COMMA int*
                                                COMMA int
                                                COMMA complex_t*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        goertzel_cmplx,              complex_t*
                                                COMMA int
                                                COMMA int*
                                                COMMA int
                                                COMMA complex_t*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        group_delay,                 double*          b
                                                COMMA double*          a
                                                COMMA int              ord
                                                COMMA int              flag
                                                COMMA double*          w
                                                COMMA int              n
                                                COMMA double*          tau);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        histogram,                   double*          x
                                                COMMA int              n
                                                COMMA int              nh
                                                COMMA double*          pedges
                                                COMMA double*          ph);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        histogram_norm,              double*          y
                                                COMMA int              n
                                                COMMA int              nh
                                                COMMA double*          x
                                                COMMA double*          w);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        idft_cmplx,                  complex_t*
                                                COMMA int
                                                COMMA complex_t*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        ifft_cmplx,                  complex_t*
                                                COMMA int
                                                COMMA fft_t*
                                                COMMA complex_t* );
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        iir,                         double           rp
                                                COMMA double           rs
                                                COMMA int              ord
                                                COMMA double           w0
                                                COMMA double           w1
                                                COMMA int              type
                                                COMMA double*          b
                                                COMMA double*          a);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        linspace,                    double
                                                COMMA double
                                                COMMA int
                                                COMMA int
                                                COMMA double*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        log_cmplx,                   complex_t*
                                                COMMA int
                                                COMMA complex_t*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        logspace,                    double
                                                COMMA double
                                                COMMA int
                                                COMMA int
                                                COMMA double*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        low2bp,                      double*          b
                                                COMMA double*          a
                                                COMMA int              ord
                                                COMMA double           w0
                                                COMMA double           wpl
                                                COMMA double           wph
                                                COMMA double*          beta
                                                COMMA double*          alpha);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        low2bs,                      double*          b
                                                COMMA double*          a
                                                COMMA int              ord
                                                COMMA double           w0
                                                COMMA double           wsl
                                                COMMA double           wsh
                                                COMMA double*          beta
                                                COMMA double*          alpha);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        low2high,                    double*          b
                                                COMMA double*          a
                                                COMMA int              ord
                                                COMMA double           w0
                                                COMMA double           w1
                                                COMMA double*          beta
                                                COMMA double*          alpha);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        low2low,                     double*          b
                                                COMMA double*          a
                                                COMMA int              ord
                                                COMMA double           w0
                                                COMMA double           w1
                                                COMMA double*          beta
                                                COMMA double*          alpha);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        matrix_eig_cmplx,            complex_t*       a
                                                COMMA int              n
                                                COMMA complex_t*       v
                                                COMMA int*             info);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        matrix_eye,                  double*          a
                                                COMMA int              n
                                                COMMA int              m);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        matrix_eye_cmplx,            complex_t*       a
                                                COMMA int              n
                                                COMMA int              m);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        matrix_mul,                  double*          a
                                                COMMA int              na
                                                COMMA int              ma
                                                COMMA double*          b
                                                COMMA int              nb
                                                COMMA int              mb
                                                COMMA double*          c);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        matrix_pinv,                 double*          a
                                                COMMA int              n
                                                COMMA int              m
                                                COMMA double*          tol
                                                COMMA double*          inv
                                                COMMA int*             info);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        matrix_print,                double*          a
                                                COMMA int              n
                                                COMMA int              m
                                                COMMA const char*      name
                                                COMMA const char*      format);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        matrix_print_cmplx,          complex_t*       a
                                                COMMA int              n
                                                COMMA int              m
                                                COMMA const char*      name
                                                COMMA const char*      format);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        matrix_svd,                  double*          a
                                                COMMA int              n
                                                COMMA int              m
                                                COMMA double*          u
                                                COMMA double*          s
                                                COMMA double*          vt
                                                COMMA int*             info);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        matrix_transpose,            double*          a
                                                COMMA int              n
                                                COMMA int              m
                                                COMMA double*          b);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        matrix_transpose_cmplx,      complex_t*       a
                                                COMMA int              n
                                                COMMA int              m
                                                COMMA complex_t*       b);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        matrix_transpose_hermite,    complex_t*       a
                                                COMMA int              n
                                                COMMA int              m
                                                COMMA complex_t*       b);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        mean,                        double*          x
                                                COMMA int              n
                                                COMMA double*          m);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        mean_cmplx,                  complex_t*       x
                                                COMMA int              n
                                                COMMA complex_t*       m);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        minmax,                      double*          x
                                                COMMA int              n
                                                COMMA double*          xmin
                                                COMMA double*          xmax);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        ones,                        double*          x
                                                COMMA int              n);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        phase_delay,                 double*          b
                                                COMMA double*          a
                                                COMMA int              ord
                                                COMMA int              flag
                                                COMMA double*          w
                                                COMMA int              n
                                                COMMA double*          tau);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        poly_z2a_cmplx,              complex_t*
                                                COMMA int
                                                COMMA int
                                                COMMA complex_t*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        polyroots,                   double*          a
                                                COMMA int              ord
                                                COMMA complex_t*       r
                                                COMMA int*             info);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        polyval,                     double*
                                                COMMA int
                                                COMMA double*
                                                COMMA int
                                                COMMA double*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        polyval_cmplx,               complex_t*
                                                COMMA int
                                                COMMA complex_t*
                                                COMMA int
                                                COMMA complex_t*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        psd_bartlett,                double*         x
                                                COMMA int             n
                                                COMMA int             nfft
                                                COMMA fft_t*          pfft
                                                COMMA double          fs
                                                COMMA int             flag 
                                                COMMA double*         ppsd 
                                                COMMA double*         pfrq);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        psd_bartlett_cmplx,          complex_t*      x
                                                COMMA int             n
                                                COMMA int             nfft
                                                COMMA fft_t*          pfft
                                                COMMA double          fs
                                                COMMA int             flag 
                                                COMMA double*         ppsd 
                                                COMMA double*         pfrq);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        psd_periodogram,             double*         x
                                                COMMA int             n
                                                COMMA int             win_type
                                                COMMA double          win_param
                                                COMMA fft_t*          pfft
                                                COMMA double          fs
                                                COMMA int             flag
                                                COMMA double*         ppsd
                                                COMMA double*         pfrq);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        psd_periodogram_cmplx,       complex_t*      x
                                                COMMA int             n
                                                COMMA int             win_type
                                                COMMA double          win_param
                                                COMMA fft_t*          pfft
                                                COMMA double          fs
                                                COMMA int             flag
                                                COMMA double*         ppsd
                                                COMMA double*         pfrq);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        psd_welch,                   double*         x
                                                COMMA int             n
                                                COMMA int             win_type
                                                COMMA double          win_param
                                                COMMA int             npsd
                                                COMMA int             noverlap
                                                COMMA fft_t*          pfft
                                                COMMA double          fs
                                                COMMA int             flag
                                                COMMA double*         ppsd
                                                COMMA double*         pfrq);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,       psd_welch_cmplx,              complex_t*      x
                                                COMMA int             n
                                                COMMA int             win_type
                                                COMMA double          win_param
                                                COMMA int             npsd
                                                COMMA int             noverlap
                                                COMMA fft_t*          pfft
                                                COMMA double          fs
                                                COMMA int             flag
                                                COMMA double*         ppsd
                                                COMMA double*         pfrq);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        randb,                       double*          x
                                                COMMA int              n
                                                COMMA random_t*        prnd);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        randb2,                      double*          x
                                                COMMA int              n
                                                COMMA random_t*        prnd);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        randi,                       int*             x
                                                COMMA int              n
                                                COMMA int              start
                                                COMMA int              stop
                                                COMMA random_t*        prnd);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        randn,                       double*          x
                                                COMMA int              n
                                                COMMA double           mu
                                                COMMA double           sigma
                                                COMMA random_t*        prnd);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        randn_cmplx,                 complex_t*       x
                                                COMMA int              n
                                                COMMA complex_t*       mu
                                                COMMA double           sigma
                                                COMMA random_t*        prnd);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        random_init,                 random_t*        prnd
                                                COMMA int              type
                                                COMMA void*            seed);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        randu,                       double*
                                                COMMA int
                                                COMMA random_t*        prnd);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        ratcompos,                   double*          b
                                                COMMA double*          a
                                                COMMA int              n
                                                COMMA double*          c
                                                COMMA double*          d
                                                COMMA int              p
                                                COMMA double*          beta
                                                COMMA double*          alpha);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        re2cmplx,                    double*
                                                COMMA int
                                                COMMA complex_t*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        readbin,                     char*           fn
                                                COMMA void**          x
                                                COMMA int*            pn
                                                COMMA int*            pm
                                                COMMA int*            dtype);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        signal_pimp,                 double*
                                                COMMA size_t
                                                COMMA double
                                                COMMA double
                                                COMMA double
                                                COMMA double
                                                COMMA double*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        signal_saw,                  double*          t
                                                COMMA size_t           n
                                                COMMA double           amp
                                                COMMA double           dt
                                                COMMA double           period
                                                COMMA double*          y);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        sin_cmplx,                   complex_t*
                                                COMMA int
                                                COMMA complex_t*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        sinc,                        double*          x
                                                COMMA int              n
                                                COMMA double           a
                                                COMMA double*          y);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        sine_int,                    double*          x
                                                COMMA int              n
                                                COMMA double*          si);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        sqrt_cmplx,                  complex_t*
                                                COMMA int
                                                COMMA complex_t*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        stat_std,                    double*          x
                                                COMMA int              n
                                                COMMA double*          s);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        stat_std_cmplx,              complex_t*       x
                                                COMMA int              n
                                                COMMA double*          s);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        sum,                         double*          x
                                                COMMA int              n
                                                COMMA double*          s);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        sum_sqr,                     double*          x
                                                COMMA int              n
                                                COMMA double*          s);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        trapint,                     double*
                                                COMMA double*
                                                COMMA int
                                                COMMA double*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        trapint_cmplx,               double*
                                                COMMA complex_t*
                                                COMMA int
                                                COMMA complex_t*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        unwrap,                      double*
                                                COMMA int
                                                COMMA double
                                                COMMA double);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        vector_dot,                  double* x
                                                COMMA double* y
                                                COMMA int     n
                                                COMMA double* p);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        verif,                       double*          x
                                                COMMA double*          y
                                                COMMA size_t           n
                                                COMMA double           eps
                                                COMMA double*          err);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        verif_data_gen,              int              len
                                                COMMA int              type
                                                COMMA char*            fn);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        verif_cmplx,                 complex_t*       x
                                                COMMA complex_t*       y
                                                COMMA size_t           n
                                                COMMA double           eps
                                                COMMA double*          err);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(void,       verif_str,                   double*          yout
                                                COMMA int              nout
                                                COMMA char*            str_msg
                                                COMMA char*            outfn
                                                COMMA char*            logfn);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(void,       verif_str_cmplx,             complex_t*          yout
                                                COMMA int              nout
                                                COMMA char*            str_msg
                                                COMMA char*            outfn
                                                COMMA char*            logfn);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        window,                      double*          w
                                                COMMA int              n
                                                COMMA int              win_type
                                                COMMA double           param);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        writebin,                    void*
                                                COMMA int
                                                COMMA int
                                                COMMA int
                                                COMMA char*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        writetxt,                    double*
                                                COMMA double*
                                                COMMA int
                                                COMMA char*     );
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        writetxt_3d,                 double*          x
                                                COMMA int              nx
                                                COMMA double*          y
                                                COMMA int              ny
                                                COMMA double*          z
                                                COMMA char*            fn);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        writetxt_3dline,             double*          x
                                                COMMA double*          y
                                                COMMA double*          z
                                                COMMA int              n
                                                COMMA char*            fn);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        writetxt_cmplx,              complex_t*       x
                                                COMMA int              n
                                                COMMA char*            fn);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        writetxt_cmplx_im,           double*
                                                COMMA complex_t*
                                                COMMA int
                                                COMMA char*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        writetxt_cmplx_re,           double*
                                                COMMA complex_t*
                                                COMMA int
                                                COMMA char*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        writetxt_int,                int*
                                                COMMA int*
                                                COMMA int
                                                COMMA char*);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        xcorr,                       double*           x
                                                COMMA int               nx
                                                COMMA double*           y
                                                COMMA int               ny
                                                COMMA int               flag
                                                COMMA int               nr
                                                COMMA double*           r
                                                COMMA double*           t);
/*----------------------------------------------------------------------------*/
DECLARE_FUNC(int,        xcorr_cmplx,                 complex_t*        x
                                                COMMA int               nx
                                                COMMA complex_t*        y
                                                COMMA int               ny
                                                COMMA int               flag
                                                COMMA int               nr
                                                COMMA complex_t*        r
                                                COMMA double*           t);
/*----------------------------------------------------------------------------*/


#ifdef __cplusplus
  }
#endif


#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup SYS_LOADING_GROUP
\fn void* dspl_load()
\brief Perform dynamic linking and load libdspl-2.0 functions.

This function attempts to link to the library `libdspl.dll` in
Windows system and the `libdspl.so` library on the Linux system.
The library is assumed to be in the same directory as the application.
user, or the path to the library is registered in the operating path variables
system.

Upon successful binding and loading of library functions, the handle is returned
libraries, as well as in the address space of the application appear
pointers to libdspl-2.0 functions.

\note
The returned handle is of type `void *`, which can be cast on Windows
to type `HINSTANCE`. In practice, this is not necessary, because this
the type is cast to `HINSTANCE` automatically if the compiler flag is set,
indicating that the application is being built on Windows.

An example of a simple program that implements dynamic binding with DSPL-2.0.

\code
#include <stdio.h>
#include <stdlib.h>
#include "dspl.h"

int main (int argc, char* argv[])
{
    void *hdspl;          // DSPL handle 
    hdspl = dspl_load (); // Dynamic linking 
    
    // Check the pointer. If `NULL`, then the link failed 
    if (! hdspl)
    {
        printf ("libdspl loading error! \n");
        return -1;
    }
    
    // The link was successful, you can call the functions of DSPL-2.0
    
    //Before correctly terminating the application, you must unlink
    //library and clear memory.
    dspl_free(hdspl);
    
    return 0;
}
\endcode

\author Bakhurin Sergey. www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup SYS_LOADING_GROUP
\fn void* dspl_load()
\brief Произвести динамическую линковку и загрузить функции libdspl-2.0.

Данная функция производит попытку связывания с библиотекой `libdspl.dll` в 
системе Windows и с библиотекой `libdspl.so` в системе Linux. 
Предполагается, что библиотека находится в одной директории с приложением
пользователя, или путь к библиотеке прописан в переменных пути операционной
системы.

При удачном связывании и загрузке функций библиотеки возвращается хэндл 
библиотеки, а также в адресном пространстве приложения появляются 
указатели на функции libdspl-2.0. 

\note
Возвращаемый хэндл имеет тип `void*`, который в ОС Windows может быть приведен
к типу `HINSTANCE`. На практике необходимости в этом, нет, потому что данный
тип приводится к `HINSTANCE` автоматически, если выставлен флаг компилятора, 
указывающий, что сборка приложения производится в ОС Windows.

Пример простейшей программы реализующей динамическое связывание с DSPL-2.0.

\code
#include <stdio.h>
#include <stdlib.h>
#include "dspl.h"

int main(int argc, char* argv[])
{
    void* hdspl;           // DSPL хэндл
    hdspl = dspl_load();   // Динамическая линковка
    
    // Проверяем указатель. Если `NULL`, то линковка прошла неудачно
    if(!hdspl)
    {
        printf("libdspl loading error!\n");
        return -1;
    }
    
    // Линковка прошла успешно можно вызывать функции DSPL-2.0
     
    // Перед корректным завершением приложения необходимо разлинковать 
    // библиотеку и очистить память.
    dspl_free(hdspl); 
   
    return 0;
}
\endcode

\author Бахурин Сергей. www.dsplib.org
***************************************************************************** */
#endif
void* dspl_load();





#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup SYS_LOADING_GROUP
\fn void dspl_free(void* handle)
\brief Cleans up the previously linked DSPL-2.0 dynamic library.

This cross-platform function clears the library `libdspl.dll` in
Windows system and from the library `libdspl.so` on the Linux system.
After cleaning the library, all functions will become unavailable.

\param [in] handle
Handle of the previously linked DSPL-2.0 library. \n
This pointer can be `NULL`, in this case no action
are being produced.

\author Bakhurin Sergey. www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup SYS_LOADING_GROUP
\fn void dspl_free(void* handle)
\brief Очищает связанную ранее динамическую библиотеку DSPL-2.0.

Данная кроссплатформенная функция производит очистку библиотеки `libdspl.dll` в 
системе Windows и с библиотеки `libdspl.so` в системе Linux. 
После очистки библиотеки все функции станут недоступны.

\param[in] handle
Хэндл прилинкованной ранее библиотеки DSPL-2.0. \n
Данный указатель может быть `NULL`, в этом случае никакие действия не 
производятся.\n\n

\author Бахурин Сергей. www.dsplib.org
**************************************************************************** */
#endif
void  dspl_free(void* handle);



#endif /* DSPL_H */

