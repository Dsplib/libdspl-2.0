clear all; close all; clc;
addpath('octave');

x = readbin('dat/real.dat');
y = mean(x);
writebin(y, 0, 'dat/mean_real.dat');
y = std(x);
writebin(y, 0, 'dat/std_real.dat');

x = readbin('dat/complex.dat');
y = mean(x);
writebin(y, 1, 'dat/mean_cmplx.dat');
y = std(x);
writebin(y, 0, 'dat/std_cmplx.dat');

