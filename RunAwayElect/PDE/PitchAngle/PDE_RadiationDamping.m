%This code presents the analytical solution of the Radiation Damping
E = 0; C = 0; R = 1;

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
f0 = @(x)( exp(-(x-0.36).^2/sigma^2) );
% f0 = @(x)( exp(-(x-0.36).^2/sigma^2)+exp(-(x+0.36).^2/sigma^2) );
phi = @(x,t)( x.*exp(-t)./sqrt(1+(exp(-2*t)-1).*x.^2) );
ExactF = @(x,t)(...
    (phi(x,t).*(1-phi(x,t).^2))./(x.*(1-x.^2)).*f0(phi(x,t)) ...
    );
source = @(x,t)(x-x);
