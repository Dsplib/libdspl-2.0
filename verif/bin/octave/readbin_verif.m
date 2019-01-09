clear all; close all; clc;
addpath('octave/functions');

x = dspl_readbin('dat/in.dat');

dspl_writebin(x, 'dat/out.dat');
