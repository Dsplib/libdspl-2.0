clear all; close all; clc;




u = - 2*pi/7;

w = [  (  cos(u) +  cos(2*u) +   cos(3*u)) / 3 - 1;
       (2*cos(u) -  cos(2*u) -   cos(3*u)) / 3;
       (  cos(u) -2*cos(2*u) +   cos(3*u)) / 3;
       (  cos(u) +  cos(2*u) - 2*cos(3*u)) / 3;
       (  sin(u) +  sin(2*u) -   sin(3*u)) / 3;
       (2*sin(u) -  sin(2*u) +   sin(3*u)) / 3;
       (  sin(u) -2*sin(2*u) -   sin(3*u)) / 3;
       (  sin(u) +  sin(2*u) + 2*sin(3*u)) / 3;]
   
   
fprintf('DFT 7 coeff\n');       
for k = 1:length(w)
  fprintf('#define DFT7_W%d      % -.24f\n', k, w(k));
end  
         
