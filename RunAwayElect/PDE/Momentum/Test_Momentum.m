%This code presents the analytical solution of the FP
% df/dt = 1/x^2 d/dx (Ca*FunCoef1*df/dx+Cf*FunCoef3*f) 
Q = 1;

C = 1;
R = 0;
E = 0;


ExactF = @(x,t)(x.^2.*(x-4).*t);
source = @(x,t)(-x.^2.*(6*R*t.*x.^3+5*C*t.*x.^2-5*E*t.*x.^2-20*R*t.*x.^2-...
    4*C*t.*x+16*E*t.*x-x.^3-24*C*t+4.*x.^2));
ExactQ = @(x,t)(2*t*x.*(x-4)+t*x.^2);


Ca = 1;
Cf = 1;


% C terms
PDE.term1.Opt = 'Grad';
PDE.term1.FunCoef = @(x)(x.^2);%FunCoef1(x);
PDE.term1.Coef =  Cf;

PDE.term2.Opt = 'Diff';
PDE.term2.FunCoef = @(x)(x.^2);%FunCoef3(x);
PDE.term2.Coef = Ca;

% E terms
PDE.term3.Opt = 'Grad';
PDE.term3.FunCoef = @(x)(x.^2);
PDE.term3.Coef = -E;

% R terms
PDE.term4.Opt = 'Grad';
PDE.term4.FunCoef = @(x)(x.^3);
PDE.term4.Coef = R;

PDE.BC.q_L = 0; PDE.BC.q_R = 1;
PDE.BC.f_L = 1; PDE.BC.f_R = 0;

