%This code presents the analytical solution of the electric field
% and collisiions

E = 4;
C = 1;
R = 0;
A = E/C;
ff0 = @(x)(x-x+1/2);
ExactF = @(x,t)(A/(2*sinh(A))*exp(A*x));
source = @(x,t)(x-x);


