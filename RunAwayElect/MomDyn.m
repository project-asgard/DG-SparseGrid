% main code for LDG 
% of Momentum Dynamics ::
% df/dt = 1/p^2*d/dp[p^2*(Ca*df/dp+Cf*f)]
%        -E/p^2*d/dp[p^2*f]
%        +R*1/p^2*d/dp*[p^3*gamma*f]
% with f(t=0)=f0(x), df/dp(p = 0)=0, f(p=pmax)=0

% Denote q = p^2*(Ca*df/dp)
% p^2 df/dt = d/dp[p^2*(Ca*df/dp+Cf*f)]
% Term 1: d/dp [Cf*p^2*f]
% Term 2: d/dp [Ca*p^2*df/dp]

clear 
close all
% clc

E = 1;
C = 1;
R = 0;

Cf = 1; Ca = 1;

% Test
PDE.term1.Opt = 'Grad';
PDE.term1.FunCoef = @(x)( (x.^2) );
PDE.term1.Coef =  Cf;

PDE.term2.Opt = 'Diff';
PDE.term2.FunCoef = @(x)( (x.^2) );
PDE.term2.Coef = Ca;

PDE_FP;

format short e
addpath(genpath(pwd))

Lev = 4;
Deg = 4;
num_plot = 3;

LInt = 0;
LEnd = 4;
Lmax = LEnd-LInt;

%% Matrix
% Term 1
Mat_Term1 = MatrixGrad(Lev,Deg,LInt,LEnd,1,PDE.term1.FunCoef);

% Term 2
Mat_Term2 = MatrixDiff(Lev,Deg,LInt,LEnd,PDE.term2.FunCoef);

MatMass = MatrixMass(Lev,Deg,LInt,LEnd,@(x)(x.^2));

% Assemble all terms
Mat = ... 
    PDE.term1.Coef * Mat_Term1 + ...
    PDE.term2.Coef * Mat_Term2 ;

%% RHS
time = 0;

rhs = ComputRHS(Lev,Deg,LInt,LEnd,source,time);

%% B.C

[x_node,Meval] = PlotDGData(Lev,Deg,LInt,LEnd);


%% Solve
DoFs = size(Mat,1);
f0 = ComputRHS(Lev,Deg,LInt,LEnd,ExactF,0);
fn = f0;
dt = ((LEnd - LInt)/2^Lev)^Deg*0.01;
for Iter = 1 : 0%1000
    
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
