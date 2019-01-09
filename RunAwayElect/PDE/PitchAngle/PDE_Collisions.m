%This code presents the analytical solution of the collision
E = 0; C = 1; R = 0;

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

% Analytical Solution
h0 = 3; h1 = 0.5; h2 = 1; h3 = 0.7; h4 = 3; h6 = 3;
ExactF = @(x,t)(h0+h1*x+h2*(3*x.^2-1)/2+h3*(5*x.^3-3*x)/2+...
    h4*(35*x.^4-30*x.^2+3)/8+h6*(231*x.^6-315*x.^4+105*x.^2-5)/16);
source = @(x,t)(x-x);



