clear all; close all; clc;
addpath('octave');

x = readbin('dat/real.dat');
writebin(x, 0, 'dat/yreal.dat');

x = readbin('dat/complex.dat');
writebin(x, 1, 'dat/ycomplex.dat');


