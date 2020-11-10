clear all; close all; clc;
addpath('octave');


fn_in = {'dat/x_fft_4.dat';
         'dat/x_fft_8.dat';
         'dat/x_fft_16.dat';
         'dat/x_fft_32.dat';
         'dat/x_fft_64.dat';
         'dat/x_fft_128.dat';
         'dat/x_fft_256.dat';
         'dat/x_fft_512.dat';
         'dat/x_fft_1024.dat';
         'dat/x_fft_2048.dat';
         'dat/x_fft_4096.dat';
         'dat/x_fft_8192.dat';
         'dat/x_fft_16384.dat';
         'dat/x_fft_32768.dat';
         'dat/x_fft_65536.dat'};

fn_out = {'dat/y_fft_4.dat';
          'dat/y_fft_8.dat';
          'dat/y_fft_16.dat';
          'dat/y_fft_32.dat';
          'dat/y_fft_64.dat';
          'dat/y_fft_128.dat';
          'dat/y_fft_256.dat';
          'dat/y_fft_512.dat';
          'dat/y_fft_1024.dat';
          'dat/y_fft_2048.dat';
          'dat/y_fft_4096.dat';
          'dat/y_fft_8192.dat';
          'dat/y_fft_16384.dat';
          'dat/y_fft_32768.dat';
          'dat/y_fft_65536.dat'};


for i = 1:length(fn_in)
    x = readbin(fn_in{i});
    y = fft(x);
    writebin(y, 1, fn_out{i});
end


