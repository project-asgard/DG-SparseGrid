%This code presents the analytical solution of the electric field
%acceleration

sigma = 0.1;
% f0 = @(x)( exp(-x.^2/sigma^2) );
f0 = @(x)(x-x+1);
phi = @(x,t)( tanh(atanh(x)-t) );
ExactF = @(x,t)( (1-phi(x,t).^2)./(1-x.^2).*f0(phi(x,t)) );
source = @(x,t)(x-x);

E = 1;
C = 0;
R = 0;