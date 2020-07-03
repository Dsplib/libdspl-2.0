# libdspl-2.0 free digital signal processing algorithm library

[![Build Status](https://travis-ci.org/Dsplib/libdspl-2.0.svg?branch=master)](https://travis-ci.org/Dsplib/libdspl-2.0)
![GitHub](https://img.shields.io/github/license/Dsplib/libdspl-2.0)

libdspl-2.0 — opensource cross-platform digital signal processing algorithm library, written in C language.
Distributed under LGPL v3 license. This allows to use this library in all applications with dynamic linking.

libdspl-2.0 includes follow algorithms sets:
* Digital spectral analysis, discrete and fast Fourier transform algorithms.
* Analog and digital IIR filters design and analysis.
* Digital FIR filters design and analysis.
* Windows function collection includes 15 different parametric and nonparametric window functions.
* Digital Hilbert transform algorithms.
* Mathematical sections includes trigonometric, hyperbolic, elliptic functions of real and complex variables.
* Pseudorandom numbers generation algorithms.
* Statistic functions.
* Linear algebra [BLAS](http://www.netlib.org/blas/) and [LAPACK](http://www.netlib.org/lapack/) packages are used under the hood.
* Digital resampling algorithms.


### Build and run libdspl-2.0
To build the DSPL-2.0 library on Windows, a special set of programs _dsplib_ _toolchain_ is provided. Dsplib toolchain includes GCC, Gnuplot, CodeBlocks IDE, file manager Far and also Unix utilities for Windows OS.


### Documentation content
* Mathematical sections:
  * [Basic math functions of the real and complex arguments.](http://en.dsplib.org/dspl/group___s_p_e_c___m_a_t_h___c_o_m_m_o_n___g_r_o_u_p.html)
  * [Trigonometric and hyperbolic of functions the real and complex arguments.](http://en.dsplib.org/dspl/group___s_p_e_c___m_a_t_h___t_r_i_g___g_r_o_u_p.html)
  * [Transcendent math functions.](http://en.dsplib.org/dspl/group___s_p_e_c___m_a_t_h___t_r_a_n_s_c_e_n_d.html)
  * [Elliptic Jacoby functions of the real and complex arguments.](http://en.dsplib.org/dspl/group___s_p_e_c___m_a_t_h___e_l_l_i_p___g_r_o_u_p.html)
  * [Pseudo-random numbers generation.](http://en.dsplib.org/dspl/group___s_p_e_c___m_a_t_h___r_a_n_d___g_e_n___g_r_o_u_p.html)
  * [Math statistic functions.](http://en.dsplib.org/dspl/group___s_p_e_c___m_a_t_h___s_t_a_t___g_r_o_u_p.html)
  * [Linear algebra and matrix operations.](http://en.dsplib.org/dspl/group___s_p_e_c___m_a_t_h___l_i_n_a_l_g___g_r_o_u_p.html)

* Digital spectral analysis:
  * [Discrete Fourier transform and fast Fourier transform algorithms.](http://en.dsplib.org/dspl/group___d_f_t___g_r_o_u_p.html)
  * [Windows function for filter design and spectrum analysis.](http://en.dsplib.org/dspl/group___w_i_n___g_r_o_u_p.html)
  * [Hilbert transform algorithms.](http://en.dsplib.org/dspl/group___h_i_l_b_e_r_t___g_r_o_u_p.html)

* Analog and digital filters design and analysis:
  * [Convolution and digital filtration.](http://en.dsplib.org/dspl/group___f_i_l_t_e_r___c_o_n_v___g_r_o_u_p.html)
  * [IIR filters design.](http://en.dsplib.org/dspl/group___i_i_r___f_i_l_t_e_r___d_e_s_i_g_n___g_r_o_u_p.html)
  * [FIR filter design.](http://en.dsplib.org/dspl/group___f_i_r___f_i_l_t_e_r___d_e_s_i_g_n___g_r_o_u_p.html)
  * [Analog and digital filter analysis.](http://en.dsplib.org/dspl/group___f_i_l_t_e_r___a_n_a_l_y_s_i_s___g_r_o_u_p.html)

* Other algorithms:
  * [Digital samplerate conversion (resampling).](http://en.dsplib.org/dspl/group___r_e_s_a_m_p_l_i_n_g___g_r_o_u_p.html)
  * [Input and output data from external files.](http://en.dsplib.org/dspl/group___i_n___o_u_t___g_r_o_u_p.html)
  * [Plotting data by Gnuplot interface.](http://en.dsplib.org/dspl/group___p_l_o_t___g_r_o_u_p.html)

* Appendix
  * [libdspl-2.0 data types.](http://en.dsplib.org/dspl/group___t_y_p_e_s___g_r_o_u_p.html)
  * [Error codes.](http://en.dsplib.org/dspl/group___e_r_r_o_r___c_o_d_e___g_r_o_u_p.html)