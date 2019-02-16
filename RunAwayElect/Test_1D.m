% test for Rhea's problem 

clear
clc
close all

format short e
addpath(genpath(pwd))

% % Test
% Eq (1).
% dn/dt + d/dx[n*V] = S 
% PDE.term1.Opt = 'Grad';

%bc = 0 :: Dirichlet Boundary
%bc = 1 :: Neumann Boundary
PDE.BC.f_L = 0; PDE.BC.f_R = 1;

Lev = 5;
Deg = 2;
num_plot = 2;

DoFs = 2^Lev*Deg;

LInt = -1;
LEnd = 1;
Lmax = LEnd-LInt;
dx = Lmax/2^Lev;


% Plotting Data
[quad_x,quad_w]=lgwt(num_plot,-1,1);
ww = repmat(quad_w,2^Lev,1)*dx/2;%2^(-Lev-1);
[x_node,Meval] = PlotDGData(Lev,Deg,LInt,LEnd,num_plot);

% Begin with Eq (1)
CFL = 0.1;
dt = CFL*(dx);
MaxT = ceil(1e-1/dt);


Eq1TestType = 1;

switch Eq1TestType
    case 1
% Test 1
ExactF = @(x,t)(t*sin(pi*x));
Source = @(x,t)(sin(pi*x)+t*sin(pi*x)+pi*t*x.*cos(pi*x));
CoeFun = @(x)(x);
BCFunc = @(x,t)(x.*t.*sin(pi*x));
    case 2
% Test 2
ExactF = @(x,t)(exp(-t)*exp(x));
Source = @(x,t)(x-x);
CoeFun = @(x)(1);
BCFunc = @(x,t)(exp(-t)*exp(x));
    case 3
% Test 3
ExactF = @(x,t)(t*sin(pi*x));
Source = @(x,t)(sin(pi*x)+pi*t*cos(pi*x));
CoeFun = @(x)(x-x+1);
BCFunc = @(x,t)(t*sin(pi*x));
end

% Term 1: d/dx[FunCoef*f]
FluxValue1 = 1; % choose upwind flux
Mat_Term1 = MatrixGrad(Lev,Deg,LInt,LEnd,FluxValue1,CoeFun, @(x)0, PDE.BC.f_L,PDE.BC.f_R);


time = 0;
n0 = zeros(DoFs,1);
n0 = ComputRHS(Lev,Deg,LInt,LEnd,ExactF,time);

for T = 1 : 1%MaxT
    rhs = ComputRHS(Lev,Deg,LInt,LEnd,Source,time);

    time = time + dt;
    
    bc = ComputeBC(Lev,Deg,LInt,LEnd,BCFunc,time,PDE.BC.f_L,PDE.BC.f_R);
    
    n1 = n0 - dt*Mat_Term1*n0 + dt*(rhs - bc);
    n0 = n1;
    
    n0_tmp = Meval*n0;
    plot(x_node,n0_tmp,'r-',x_node,ExactF(x_node,time),'b--','LineWidth',2)
    title(['time = ', num2str(time)])
    
    % check about Mass and Engergy
    Mas_tol(T) = sum(ww.*n0_tmp);
    Eng_tol(T) = sum(ww.*n0_tmp.^2);
    pause(0.1)
end

figure;plot(Mas_tol,'r-','LineWidth',2);hold on;plot(Eng_tol,'b-','LineWidth',2)
title('Total Mass');

['Begin with Eq 2']
% Begin with Eq (2)
% Eq(2): dV/dt + V*dV/dx - Gam*d^2V/dx^2 = f
% Let Q = dV/dx
% => dV/dt + V*dV/dx - Gam*dQ/dx = f
%    Q = dV/dx

CFL = .1;%.1;
dt = CFL*(dx);%^Deg;
MaxT = ceil(1e1/dt);

gamma = 0.001;%-1e-0; % gamma = -1 is correct, but fail on other cases

% Test
ExactF = @(x,t)(1-tanh(1/(2*gamma)*(x-t)));
ExactQ = @(x,t)( -(1/2)*(1-tanh((1/2)*(x-t)/gamma).^2)/gamma );
Source = @(x,t)(x-x);
BCFunc = @(x,t)(1-tanh(1/(2*gamma)*(x-t)));
BCBurgers = @(x,t)( (1-tanh(1/(2*gamma)*(x-t))).^2);

% Test 
% ExactF = @(x,t)(sin(pi*x));
% ExactQ = @(x,t)(pi*cos(pi*x));
% Source = @(x,t)(x-x);
% BCFunc = @(x,t)(sin(pi*x));
% BCBurgers = @(x,t)( (sin(pi*x)).^2);

% % % Test
% % beta = 0.5;
% % ExactF = @(x,t)(beta*sin(x));
% % ExactQ = @(x,t)(beta*cos(x));
% % Source = @(x,t)(sin(x).*cos(x));
% % BCFunc = @(x,t)(x-x);

% % Burgers Test 2
% ExactF = @(x,t)(1+sin(pi*x));
% ExactQ = @(x,t)(x-x);
% Source = @(x,t)(x-x);
% BCFunc = ExactQ;
% BCBurgers = @(x,t)((1+sin(pi*x)).^2/2);

f_bcL = 0; f_bcR = 0;
q_bcL = 0; q_bcR = 0;

Mat1 =  MatrixGradU(Lev,Deg,LInt,LEnd,1,@(x)1,@(x)0,f_bcL,f_bcR); % equation for q
Mat2 = -MatrixGradQ(Lev,Deg,LInt,LEnd,-1,@(x)1,@(x)0,q_bcL,q_bcR); % equation for f

time = 0;
% Initial Condition
n0 = ComputRHS(Lev,Deg,LInt,LEnd,ExactF,time);
nn0 = n0;
FluxValue = 1;
figure;

BurbcL = 0; BurbcR = 1;

for T = 1 : 100%MaxT
    time = time + dt;
    rhs = ComputRHS(Lev,Deg,LInt,LEnd,Source,time);
    % B.C for Diffusion Part
    bc = ComputeBC(Lev,Deg,LInt,LEnd,BCFunc,time,f_bcL,f_bcR,FluxValue);

% %     % Burgers Part
% %     [uL,uR] = BurSign(LInt,LEnd,Lev,Deg,n0);
% %     if uL >=0
% %         BurbcL = 0;
% %     else
% %         BurbcL = 1;
% %     end
% %     if uR <=0
% %         BurbcR = 0;
% %     else
% %         BurbcR = 1;
% %     end
    [FF,MaxC] = VectorBurgers(Lev,Deg,LInt,LEnd, 1,n0,BurbcL,BurbcR);% equation for q
    % B.C for Burgers Part
    bc1 = ComputeBCBurgers(Lev,Deg,LInt,LEnd,BCFunc,time,BurbcL,BurbcR,MaxC);
    
    % Explicit Advance
%     n1 = n0  + dt*(rhs) - dt*FF - dt*bc1 - gamma*dt*Mat2*bc -gamma* dt*(Mat2*Mat1)*n0 ;
%     n0 = n1;
    
    % Implicit Advance
    Mat = speye(DoFs,DoFs) + gamma*dt*(Mat2*Mat1);
    Mat = inv(Mat);
    n1 = Mat*n0 + Mat*(dt*rhs - gamma*dt*Mat2*bc )- dt*Mat*(FF +bc1);
    n0 = n1;
    
    n0_tmp = Meval*n0;
    
    % check about Mass and Engergy
    Mas_tol(T) = sum(ww.*n0_tmp);
    Eng_tol(T) = sum(ww.*n0_tmp.^2);
    pause(0.1)
%     
    subplot(1,2,1)
    plot(x_node,ExactF(x_node,0),'k-',x_node,Meval*n0,'r-',x_node,ExactF(x_node,time),'b--','LineWidth',2)
    title(['time = ', num2str(time)])
    subplot(1,2,2)
    plot(x_node,Meval*(Mat1*n0+bc),'r-',x_node,ExactQ(x_node,time),'b--','LineWidth',2)
    title(['time = ', num2str(time)])
    
    sol(:,T) = Meval*n0;
    pause(0.1)
end

figure;
plot(Mas_tol,'ro');hold on; plot(Eng_tol,'b+')
 