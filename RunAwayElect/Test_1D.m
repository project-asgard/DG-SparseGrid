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


LInt = 0;
LEnd = .1;
Lmax = LEnd-LInt;
dx = Lmax/2^Lev;

CFL = 0.01;
dt = CFL*(dx);
MaxT = ceil(1e-1/dt);

% Test 1
ExactF = @(x,t)(t*sin(pi*x));
Source = @(x,t)(sin(pi*x)+t*sin(pi*x)+pi*t*x.*cos(pi*x));
CoeFun = @(x)(x);
BCFunc = @(x,t)(x.*t.*sin(pi*x));

% Test 2
% ExactF = @(x,t)(t*sin(pi*x));
% Source = @(x,t)(sin(pi*x)+pi*t*cos(pi*x));
% CoeFun = @(x)(x-x+1);
% BCFunc = @(x,t)(t*sin(pi*x));

% Term 1: d/dx[FunCoef*f]
FluxValue1 = 1; % choose upwind flux
Mat_Term1 = MatrixGrad(Lev,Deg,LInt,LEnd,FluxValue1,CoeFun, @(x)0, PDE.BC.f_L,PDE.BC.f_R);

DoFs = 2^Lev*Deg;
n0 = zeros(DoFs,1);

time = 0;

[quad_x,quad_w]=lgwt(num_plot,-1,1);
ww = repmat(quad_w,2^Lev,1)*2^(-Lev-1);

[x_node,Meval] = PlotDGData(Lev,Deg,LInt,LEnd,num_plot);


n0 = ComputRHS(Lev,Deg,LInt,LEnd,ExactF,time);

for T = 1 : MaxT
    rhs = ComputRHS(Lev,Deg,LInt,LEnd,Source,time);

    time = time + dt;
    
    bc = ComputeBC(Lev,Deg,LInt,LEnd,BCFunc,time,PDE.BC.f_L,PDE.BC.f_R);
    
    n1 = n0 - dt*Mat_Term1*n0 + dt*(rhs - bc);
    n0 = n1;
    
    n0_tmp = Meval*n0;
    plot(x_node,n0_tmp,'r-',x_node,ExactF(x_node,time),'b--')
    title(['time = ', num2str(time)])
    
    % check about Mass 
    Mass_tol(T) = sum(ww.*n0_tmp);
    
    pause(0.1)
end

figure;plot(Mass_tol);title('Total Mass');

return
% Eq(2): dV/dt + V*dV/dx - Gam*d^2V/dx^2 = f
% Let Q = dV/dx
% => dV/dt + V*dV/dx - Gam*dQ/dx = f
%    Q = dV/dx

CFL = 0.01;
dt = CFL*(dx)^3;
MaxT = ceil(1e-1/dt);

Gam = 1;
f_bcL = 0; f_bcR = 0;
q_bcL = 1; q_bcR = 1;

Mat1 = MatrixGrad(Lev,Deg,LInt,LEnd,-1,@(x)1,@(x)0,f_bcL,f_bcR); % equation for q
Mat2 = MatrixGrad(Lev,Deg,LInt,LEnd, 1,@(x)1,@(x)0,q_bcL,q_bcR); % equation for f

Mat2 =  Mat2*Mat1;

% Mat1 = MatrixGrad(Lev,Deg,LInt,LEnd, 1,@(x)(1),@(x)0,f_bcL,f_bcR);% equation for q

% Test
ExactF = @(x,t)( t*sin(pi*x) );
Source = @(x,t)( sin(pi*x)-pi^2*t.*sin(pi*x)+t*sin(pi*x).*t.*pi.*cos(pi*x) );
% CoeFun = @(x)(x);
BCFunc = @(x,t)( t.*sin(pi*x) );

[x_node,Meval] = PlotDGData(Lev,Deg,LInt,LEnd,num_plot);

n0 = ComputRHS(Lev,Deg,LInt,LEnd,ExactF,time);

for T = 1 : MaxT
    rhs = ComputRHS(Lev,Deg,LInt,LEnd,Source,time);

    time = time + dt;
    
    bc = ComputeBC(Lev,Deg,LInt,LEnd,BCFunc,time,f_bcL,f_bcR);
    
    Mat1 = MatrixGradTD(Lev,Deg,LInt,LEnd, 1,n0,f_bcL,f_bcR);% equation for q
    
    bc1 = ComputeBC(Lev,Deg,LInt,LEnd,@(x,t)(t*sin(pi*x).*t.*pi.*cos(pi*x)),time,f_bcL,f_bcR);
    
    n1 = n0 - dt*(Mat1+Mat2)*n0 + dt*(rhs - bc - bc1);
    n0 = n1;
    
    
    plot(x_node,Meval*n0,'r-',x_node,ExactF(x_node,time),'b--')
    title(['time = ', num2str(time)])
    
    pause(0.1)
end

