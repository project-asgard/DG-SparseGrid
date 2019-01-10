% main code for Localized-DG 
% of Full Pitch Angle Dynamics ::
% df/dt =  E*{-d/dx[(1-x^2)f]}+C*d/dx[(1-x^2)df/dx]-R*d/dx[x(1-x^2)f]
% with f(t=0)=f0(x), f(x=+,-1)=0.

% clc
clear 
close all

format short e
addpath(genpath(pwd))

E = 1; C = 1; R = 3;

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




% PDE_ElectricFieldAcceleration;
% PDE_Collisions;
% PDE_RadiationDamping;
% PDE_ElectricFieldCollisions;
PDE_FullPitchAngle;

Lev = 5;
Deg = 3;
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
f0 = ComputRHS(Lev,Deg,LInt,LEnd,ExactF,0);%zeros(DoFs,1);
dt = ((LEnd - LInt)/2^Lev)^Deg*0.01;
for Iter = 1 : 1000
    
    time = dt*Iter;
    rhs = ComputRHS(Lev,Deg,LInt,LEnd,source,time);
    fn = f0 + dt * Mat*f0 + dt * rhs;
    f0 = fn;
    
    plot(x_node,Meval*fn,'r-o',x_node,ExactF(x_node,time),'b--','LineWidth',2)
    title(num2str(time))
    pause(0.01)
    
end


figure;
plot(x_node,Meval*fn,'r-o',x_node,ExactF(x_node,time),'b--','LineWidth',2)
legend('Numerical Solution','Exact Solution')

val = Meval*fn - ExactF(x_node,time);

max(abs(val))