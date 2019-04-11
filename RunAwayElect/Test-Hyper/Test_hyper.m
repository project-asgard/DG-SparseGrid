% main code for Localized-DG 
% of Full Pitch Angle Dynamics ::
% df/dt =  E*{-d/dx[(1-x^2)f]}+C*d/dx[(1-x^2)df/dx]-R*d/dx[x(1-x^2)f]
% with f(t=0)=f0(x), f(x=+,-1)=0.

% clc
clear 
close all

format short e
addpath(genpath(pwd))


PDE_ElectricFieldAcceleration;
% PDE_Collisions;
% PDE_RadiationDamping;
% PDE_ElectricFieldCollisions;
% PDE_FullPitchAngle;

Lev = 4;
Deg = 2;
num_plot = Deg;

LInt = -1;
LEnd = 1;
Lmax = LEnd-LInt;

%% Matrix
% Term 1
Mat = MatrixGrad2(Lev,Deg,LInt,LEnd,0,@(x)1);



%% RHS
time = 0;
source = @(x,t)(cos(t).*cos(2*pi*x)+2*pi*sin(2*pi*x).*sin(t));
rhs = ComputRHS(Lev,Deg,LInt,LEnd,source,time);

%% B.C

[x_node,Meval] = PlotDGData(Lev,Deg,LInt,LEnd,num_plot);


%% Solve
DoFs = size(Mat,1);
ExactF = @(x,t)(sin(t).*cos(2*pi*x))
f0 = ComputRHS(Lev,Deg,LInt,LEnd,ExactF,0);
% for damping
% CFL = 1;
% 2nd order
% CFL = 0.02;
CFL = 0.01;

dx = ((LEnd - LInt)/2^Lev);
% dt = ((LEnd - LInt)/2^Lev)^(2*(Deg)/3)*CFL;
dt = CFL*dx;%*2^(2-Lev);
MaxIter = ceil(1/dt);

[quad_x,quad_w]=lgwt(Deg,-1,1);
quad_w = 2^(-Lev)/2*quad_w;
ww = repmat(quad_w,2^Lev,1);

for Iter = 1 : MaxIter
    
%     time = dt*Iter;
%     rhs = ComputRHS(Lev,Deg,LInt,LEnd,source,time);
%     fn = f0 + dt * Mat*f0 + dt * rhs;
%     f0 = fn;
    
    % 3rd-RK
    time = dt*Iter;
    rhs0 = ComputRHS(Lev,Deg,LInt,LEnd,source,time-dt);
    rhs =  ComputRHS(Lev,Deg,LInt,LEnd,source,time);
    rhs2 = ComputRHS(Lev,Deg,LInt,LEnd,source,time-dt/2);
    f1 = f0 + dt*(  Mat*f0 +rhs0 );
    f2 = 3/4*f0+1/4*f1+1/4*dt*(Mat*f1+rhs);
    fn = 1/3*f0+2/3*f2+2/3*dt*(Mat*f2+rhs2);
    f0 = fn;
    
%     qn = Mat1*fn;
    
    val = Meval*fn;
    Exval = ExactF(x_node,time);
    plot(x_node,val,'r-o',x_node,Exval,'b--',...
         ...x_node,Meval*qn,'g-o',...
        'LineWidth',2)
    title(num2str(time))
    pause(0.01)
    max(abs(Exval-val))
%     if time>2%abs(time-0.5)<dt || abs(time-1)<dt || abs(time-2)<dt || abs(time-3)<dt
% % %          pause
% 1111
%     end
    TolTime(Iter) = time;
    TolMass(Iter) = sum(ww.*val);
    TolExMass(Iter) = sum(ww.*Exval);
    TolEng(Iter) = sqrt(sum(ww.*val.^2));
    TolExEng(Iter) = sqrt(sum(ww.*Exval.^2));
end


figure;
plot(x_node,Meval*fn,'r-o',x_node,ExactF(x_node,time),'b--','LineWidth',2)
legend('Numerical Solution','Exact Solution')

val = Meval*fn - ExactF(x_node,time);

ErrMax = max(abs(val));

[quad_x,quad_w]=lgwt(Deg,-1,1);
quad_w = 2^(-Lev)/2*quad_w;
ww = repmat(quad_w,2^Lev,1);
ErrL2 = sqrt(sum(ww.*val.^2));

[ErrMax ErrL2]
