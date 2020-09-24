clear all; close all; clc;
addpath('octave');

x = readbin('dat/x_fft_double_radix2.dat');
y = fft(x);
writebin(y, 1, 'dat/y_fft_double_radix2.dat');


x = readbin('dat/x_fft_complex_radix2.dat');
y = fft(x);
writebin(y, 1, 'dat/y_fft_complex_radix2.dat');


x = readbin('dat/x_fft_double_common.dat');
y = fft(x);
writebin(y, 1, 'dat/y_fft_double_common.dat');


x = readbin('dat/x_fft_complex_common.dat');
y = fft(x);
writebin(y, 1, 'dat/y_fft_complex_common.dat');