%**************************************************************************
% This MATLAB-Interface
% is an example to solve Vlasov Equations with Poisson/Maxwell Equations
% Vlasov Equation with Poisson solver:
% (1). f_t + v*grad_x f + E*grad_v f=0;
% (2). -Del (phi) = 1 - int_v f dv and E= (phi)'
% Vlasov Equation with Maxwell solver:
% (1). f_t + v*grad_x f + (E+B)*grad_v f=0;
% (2). Maxwell Eq.
%**************************************************************************
clc
clear
close all
format short e

addpath(genpath(pwd))

%=============================================================
%% Step 1. Setting Initial Data and Parameters
% Lev: Level of the Mesh
%   Lev, Lev_x, and Lev_v
% k: Degree of Polynomial
% dim: dimensionality
%   dim, dim_x, and dim_v
% Lmax: x in [0,Lmax]
% Vmax: v in [-Vmax,Vmax]
% fval: initial condition for f
% rho: initial value of int_v fdv
% Maxwell coefficients: nu and eps
%=============================================================

% Ending time and test pde
T = 5;
pde = Vlasov4;
Vmax = pde.Vmax;
Lmax = pde.Lmax;

% Level information
Lev = 2;
Lev_x = Lev;Lev_v = Lev;

% Polynomial Degree
k = 1;

% Dimensionality
dim = 2;
dim_x = 1;dim_v = 1;

% Time step
dt = Lmax/2^Lev_x/Vmax/(2*k+1)/5;

%********************************
% Generate the initial condition
% Input: Lev_x,Lev_v,k,Lmax,Vmax,pde
% Output: fval--intial condition f(x,v,t=0)
%                rho--intial condition rho(x,t=0)
%                Eng--computed energy
%********************************
[fval,rho,Eng] = Intial_Con(Lev_x,Lev_v,k,Lmax,Vmax,pde);

%=============================================================
%% Step 2. Generate time-independent coefficient matrices
% Vlasolv Solver:
%   Operators:  vMassV: int_v v*l_i(v)*l_j(v)dv
%               GradV: int_v (l_i(v))'*l_j(v)dv
%               GradX: int_x (m_i(x))'*m_j(x)dx
% Poisson Solver:
%   Operators: DelaX: int_x (m_i(x))''*m_j(x)dx
% Maxwell Solver:
%   Operators: CurlCurlX: int_x curl(m_i(x))*m_j(x)dx
% Input: Lev_x, Lev_v, k, dim, Lmax, Vmax
% Output: 2D Matrices--vMassV,GradV,GradX,DeltaX, and CurlCurl(missing)
%=============================================================
[vMassV,GradV,GradX,DeltaX] = matrix_coeff_TI(Lev_x,Lev_v,k,Lmax,Vmax);


%=============================================================
%% Step 3. Generate Sparse Grids/Hash Table
% Input:    Lev-Level of the Mesh
%           k  -Degree of Polynomial
%           dim-dimensionality
% Output: Hash Table
%=============================================================
% SparseGrid;


%=============================================================
%% Step 4. Generate time-independent global Matrix
% Compute the global matrix for spacial variables "X" by
%
% Poisson Solver: A_Poisson (Hash, Dim_x,k,Lev_x,DeltaX)
% or Maxwell Solver: A_Maxwell (Hash,Dim_x,k,Lev_x,nu,eps,CurlCurlX)
% Input: Hash, Dim_x,k,Lev_x,DeltaX,or nu, eps,CurlCurlX
% Output: A_Poisson or A_Maxwell
%=============================================================
% return


%=============================================================
%% Step 5. Time Loop
%	Step 5.1 Vlasov Equation
%       Generate time dependent coefficient matrix
%       Generate global matrix A_Vlasov(Hash,coef_mat,Dim)
%       Apply A_Vlasov->f by RK
%	Step 5.2 Poisson Equation
%       Calculate rho(f)
%       Solve Poisson Equation sol_Poisson by A_Poisson(rho)
%       Compute E=(sol_Poisson)'
%or Step 5.2 Maxwell Equation
%       Calculate J(f)
%       Solve Maxwell Equation sol_Maxwell by A_Maxwell(J)
% Note:
% 3-rd order Runge-Kutta Methods
%   f1=f0+dt*(A*f0)
%   f2=3/4*f0+1/4*f1+1/4*dt*(A*f1)
%	f =1/3*f0+2/3*f2+2/3*dt*(A*f2)
% capability: vary the time-integration order of RK methods
%=============================================================
% Preparing the Plotting Data
% Plotting Data
[Meval_v,v_node,Meval_x,x_node]=matrix_plot(Lev_x,Lev_v,k,Lmax,Vmax);
[xx,vv]=meshgrid(x_node,v_node);

count=1;
for step_num = 1:floor(T/dt)
    %     time = dt*step_num;
    
    %=============================================================
    %% Step 5.3. Generate time-dependent coefficient matrix
    % Vlasolv Solver:
    %   Operators:  EMassX: int_x E(x)*m_i(x)*m_j(x)dx
    % Input: Lev, k, dim, Lmax, Vmax
    % Output: EMassX
    % Note: E is solved by Poisson or Maxwell's equation
    %=============================================================
    % Solve E from rho
    EE = Poisson_solver(Lev_x,k,Lmax,rho,DeltaX,Eng);
    % Construct EMassX
    EMassX = matrix_coeff_TD(Lev_x,k,Lmax,EE);
    
    %====================================
    % RK Time Stepping Method
    %====================================
    f_1 = fval+dt*(Multi_2D(vMassV,GradX,fval)+Multi_2D(GradV,EMassX,fval));
    f_2 = 3/4*fval+1/4*f_1+1/4*dt*(Multi_2D(vMassV,GradX,f_1)+Multi_2D(GradV,EMassX,f_1));
    fval = 1/3*fval+2/3*f_2+2/3*dt*(Multi_2D(vMassV,GradX,f_2)+Multi_2D(GradV,EMassX,f_2));
    
    % Compute rho from f(x,v,t=time)
    [rho,Eng]=Comput_rho(Lev_x,Lev_v,k,Lmax,Vmax,fval);
    
    % Check conservation properties
    
    % Plot solution
    time(count)=step_num*dt;
    rho_time(count)=abs(rho(1));
    
    count=count+1;
    
    
    % Plotting Numerical Solution
    figure(1000)
    set(gcf, 'Position', [100, 100, 1200, 900])
    subplot(1,2,1)
    tmp=Multi_2D(Meval_v,Meval_x,fval);
    
    mesh(xx,vv,reshape(tmp,k*2^Lev_x,k*2^Lev_v)','FaceColor','interp','EdgeColor','interp');
    axis([0 Lmax -Vmax Vmax])
    view(0,90)
    colorbar
    
    title(['Time at ',num2str(step_num*dt)])
    subplot(1,2,2);hold on;
    semilogy(time,rho_time,'r--');
    title([' step num = ',num2str(step_num)])
    
    pause(0.01)
    
end