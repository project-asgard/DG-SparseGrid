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

PDE_FP2;

format short e
addpath(genpath(pwd))

Lev = 4;
Deg = 3;
num_plot = 3;

LInt = 0;
LEnd = 4;
Lmax = LEnd-LInt;

%% Matrix
% Term 1: d/dx[FunCoef*f]
Mat_Term1 = MatrixGrad(Lev,Deg,LInt,LEnd,0,PDE.term1.FunCoef, @(x)0, PDE.BC.f_L,PDE.BC.f_R);

% Term 2: d/dx[FunCoef*df/dx]
[Mat_Term2,Mat1,Mat2] = MatrixDiff_Momentum(Lev,Deg,LInt,LEnd,PDE.term2.FunCoef,PDE.BC.q_L ,PDE.BC.q_R ,PDE.BC.f_L ,PDE.BC.f_R);

MatMass = MatrixMass(Lev,Deg,LInt,LEnd,@(x)(x.^2));
MatMass = inv(MatMass);

% Assemble all terms
Mat = ... 
    PDE.term1.Coef * Mat_Term1 + ...
    PDE.term2.Coef * Mat_Term2 ;

%% RHS
time = 0;

rhs = ComputRHS(Lev,Deg,LInt,LEnd,source,time);

%% B.C

[x_node,Meval] = PlotDGData(Lev,Deg,LInt,LEnd);

qq = @(x)(-2*x.*exp(-x.^2));

%% Solve
DoFs = size(Mat,1);
f0 = ComputRHS(Lev,Deg,LInt,LEnd,ExactF,0);
q0 = ComputRHS(Lev,Deg,LInt,LEnd,@(x,t)(-2*x.*exp(-x.^2)),0);
fn = f0;
dt = ((LEnd - LInt)/2^Lev)^(Deg/3)*0.001;
for Iter = 1 : 1000
    
    time = dt*Iter;
    rhs0 = ComputRHS(Lev,Deg,LInt,LEnd,source,time-dt);
    rhs = ComputRHS(Lev,Deg,LInt,LEnd,source,time);
    rhs2 = ComputRHS(Lev,Deg,LInt,LEnd,source,time+dt);
    
%     % Euler
%     fn = f0 + dt * MatMass * Mat*f0 + dt * MatMass *rhs;
%     f0 = fn;
    
    % 3rd-RK
    f1 = f0 + dt* MatMass *(  Mat*f0 +rhs0 );
    f2 = 3/4*f0+1/4*f1+1/4*dt*MatMass *(Mat*f1+rhs);
    fn = 1/3*f0+2/3*f2+2/3*dt*MatMass *(Mat*f2+rhs2);
    f0 = fn;
    
    qn = Mat1*fn;
    
    plot(x_node,Meval*fn,'r-o',x_node,ExactF(x_node,time),'b--',...
         x_node,Meval*qn,'g-<',x_node,qq(x_node),'b--',x_node,Meval*(Mat2*fn),'b-<',...
         'LineWidth',2)
    title(num2str(time))
    pause(0.01)
    
end


figure;
plot(x_node,Meval*fn,'r-o',x_node,ExactF(x_node,time),'b--','LineWidth',2)
legend('Numerical Solution','Exact Solution')

val = Meval*fn - ExactF(x_node,time);

max(abs(val))
