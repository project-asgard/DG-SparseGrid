%This code presents the analytical solution of the electric field
%acceleration
E = 1; C = 0; R = 0;

PDE.term1.Opt = 'Grad';
PDE.term1.FunCoef = @(x)( (1-x.^2) );
PDE.term1.Coef =  -E;

PDE.term2.Opt = 'Diff';
PDE.term2.FunCoef = @(x)( (1-x.^2) );
PDE.term2.Coef = C;

PDE.term3.Opt = 'Grad';
PDE.term3.FunCoef = @(x)( -(x.*(1-x.^2)) );
PDE.term3.Coef =  R;

PDE.BC.q_L = 0; PDE.BC.q_R = 0;
PDE.BC.f_L = 1;  PDE.BC.f_R = 1;

sigma = 0.1;
% f0 = @(x)( exp(-x.^2/sigma^2) );
f0 = @(x)(x-x+1);
phi = @(x,t)( tanh(atanh(x)-t) );
ExactF = @(x,t)( (1-phi(x,t).^2)./(1-x.^2).*f0(phi(x,t)) );
source = @(x,t)(x-x);

