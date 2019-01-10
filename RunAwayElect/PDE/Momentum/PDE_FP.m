%This code presents the analytical solution of the FP
Q = 1;
ExactF = @(x,t)(Q*exp(-x.^2));

source = @(x,t)(-2*x.*exp(-x.^2).*(2*x.^3-x.^2-3*x+1));
psi = @(x,t)(1./x.^2.*(erf(x)-2*x/sqrt(pi).*exp(-x.^2)));

Ca = 1;
Cf = 1;

% % Test
PDE.term1.Opt = 'Grad';
PDE.term1.FunCoef = @(x)( (x.^2) );
PDE.term1.Coef =  Cf;

PDE.term2.Opt = 'Diff';
PDE.term2.FunCoef = @(x)( (x.^2) );
PDE.term2.Coef = Ca;
