%This code presents the analytical solution of the electric field
% and collisiions
E = 4; C = 1; R = 0;

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

A = E/C;
ff0 = @(x)(x-x+1/2);
ExactF = @(x,t)(A/(2*sinh(A))*exp(A*x));
source = @(x,t)(x-x);


