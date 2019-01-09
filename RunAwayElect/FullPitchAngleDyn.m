% main code for Localized-DG 
% of Full Pitch Angle Dynamics ::
% df/dt =  E*{-d/dx[(1-x^2)f]}+C*d/dx[(1-x^2)df/dx]-R*d/dx[x(1-x^2)f]
% with f(t=0)=f0(x), f(x=+,-1)=0.

% clc
clear 
close all

format short e
addpath(genpath(pwd))

E = 1; C = 1; R = 1;

% Test 1
PDE.term1.Opt = 'Grad';
PDE.term1.FunCoef = @(x)( (1-x.^2) );
PDE.term1.Coef =  -E;

PDE.term2.Opt = 'Diff';
PDE.term2.FunCoef = @(x)( (1-x.^2) );
PDE.term2.Coef = C;

PDE.term3.Opt = 'Grad';
PDE.term3.FunCoef = @(x)( -(x.*(1-x.^2)) );
PDE.term3.Coef =  R;

% % Exact Solution
% sigma = 0.1;
% f0 = @(x)( exp(-x.^2/sigma^2) );
% % f0 = @(x)(x-x+1);
% phi = @(x,t)( tanh(atanh(x)-t) );
% ExacF = @(x,t)(...
%     (1-phi(x,t).^2)./(1-x.^2).*f0(phi(x,t)) ...
%     );

ExactF = @(x,t)(sin(pi*x)*t);
source = @(x,t)( sin(pi*x)+...
    +E*t*(-2*x.*sin(pi*x)+(-x.^2+1).*cos(pi*x)*pi)+...
    -C*t*(-2*x.*cos(pi*x)*pi-(-x.^2+1).*sin(pi*x)*pi^2)+...
    +R*t*( (-x.^2+1).*sin(pi*x)-2*x.^2.*sin(pi*x)+x.*(-x.^2+1).*cos(pi*x)*pi) );

Lev = 5;
Deg = 2;
num_plot = 3;

LInt = -1;
LEnd = 1;
Lmax = LEnd-LInt;

%% Matrix
% Term 1
Mat_Term1 = MatrixGrad(Lev,Deg,LInt,LEnd,1,PDE.term1.FunCoef);

% Term 2
Mat_Term2 = MatrixDiff(Lev,Deg,LInt,LEnd,PDE.term2.FunCoef);

% Term 3
Mat_Term3 = MatrixGrad(Lev,Deg,LInt,LEnd,1,PDE.term3.FunCoef);

% Assemble all terms
Mat = ... 
    PDE.term1.Coef * Mat_Term1 + ...
    PDE.term2.Coef * Mat_Term2 + ...
    PDE.term3.Coef * Mat_Term3 ;

%% RHS
time = 0;

rhs = ComputRHS(Lev,Deg,LInt,LEnd,source,time);

%% B.C

[x_node,Meval] = PlotDGData(Lev,Deg,LInt,LEnd);


%% Solve
DoFs = size(Mat,1);
f0 = zeros(DoFs,1);
dt = ((LEnd - LInt)/2^Lev)^Deg*0.01;
for Iter = 1 : 1000
    
    time = dt*Iter;
    rhs = ComputRHS(Lev,Deg,LInt,LEnd,source,time);
    fn = f0 + dt * Mat*f0 + dt * rhs;
    f0 = fn;
    
    plot(x_node,Meval*fn,'r-o',x_node,ExactF(x_node,time),'b--','LineWidth',2)
%     pause
end


% figure
% plot(x_node,Meval*rhs,'b--', x_node,source(x_node,time),'r--','LineWidth',2)

figure;
% checked of projection
plot(x_node,Meval*fn,'r-o',x_node,ExactF(x_node,time),'b--','LineWidth',2)

val = Meval*fn - ExactF(x_node,time);

max(abs(val))