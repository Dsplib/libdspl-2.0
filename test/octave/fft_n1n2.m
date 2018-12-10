clear all;
close all;
clc;

N  = 16;
N1 = 4;
N2 = 4;

% входной сигнал это вектор столбец размерности [N x 1]
x = (0:N-1)';

%x = [4;6;8;10];

for n1 = 0 : N1-1
  for k2 = 0 : N2-1
    W(n1 + 1, k2 + 1) = exp(-2i * pi * n1 * k2 / N);
  end
end


A = reshape(x, N2, N1);

B = A.';
D = fft(B);
F = D.*W;
G = F.';
H = fft(G);
P = H.';

y = [reshape(P,N,1), fft(x)]

