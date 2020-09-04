clear all; close all; clc;
addpath('octave');

x = readbin('dat/x.dat');
if(isreal(x))
    writebin(x, 0, 'dat/y.dat');
else
    writebin(x, 1, 'dat/y.dat');
end
