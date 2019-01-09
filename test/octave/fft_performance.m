clear all; close all; clc;

n = 8388608;
m = 4;

x0 = (0:n-1)+1i*(0:n-1);

while(n > 4)
  x = x0(1:n);
  fprintf('n = %12d      ', n);
  t0 = tic();  
  for k = 1:m
    X = fft(x);
  end
  dt = toc(t0) / m;
  fprintf('%12.6f\n', dt*1000);
  n = n/2;
  m = m*1.5;

end 