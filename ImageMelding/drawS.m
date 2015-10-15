function [T] =  drawS( L, gamma)
%L = 300;
%gamma = 2.1;
T1 = [1:L/2] / (L/2);
T1 = T1 .^ gamma;
%T2= T1;
T2 = 1- flipdim(T1,2);
T = double([T2 (1+T1)] / 2);
%plot(T);




