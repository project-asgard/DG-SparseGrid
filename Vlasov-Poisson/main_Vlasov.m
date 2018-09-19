%**************************************************************************
% This MATLAB-Interface for DG-SG Vlasov
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
%% Step 1. Setting Parameters
% Lev: Level of the Mesh
%   Lev, LevX, and LevV
% Deg: Degree of Polynomial
% Dim: dimensionality
%   Dim, Dim_x, and Dim_v
% Lmax: x in [0,Lmax]
% Vmax: v in [-Vmax,Vmax]
% fval: initial condition for f
% rho: initial value of int_v fdv
% Maxwell coefficients: nu and eps
%=============================================================
% Test PDE and Ending Time

pde = Vlasov4;
Vmax = pde.Vmax;
Lmax = pde.Lmax;
TEND = 10;

% Level information
Lev = 3;
LevX = Lev;LevV = Lev;

% Polynomial Degree
Deg = 2;% Deg = 2 Means Linear Element

% Dimensionality
Dim = 2;
DimX = 1;DimV = 1;

% Time step
dt = Lmax/2^LevX/Vmax/(2*Deg+1);

%*************************************************
%% Step 1.1. Set Up Matrices for Multi-wavelet
% Input:  Deg and Lev
% Output: Convert Matrix FMWT_COMP
%*************************************************
FMWT_COMP_x = OperatorTwoScale(Deg,2^LevX);
FMWT_COMP_v = OperatorTwoScale(Deg,2^LevV);

%*************************************************
%% Step 1.2. Initial Data according to Hash
% Generate the initial condition
% Input: LevX,LevV,Deg,Lmax,Vmax,pde
% Output: fval (fv and fx)--intial condition f(x,v,t=0)
%         rho--intial condition rho(x,t=0)
%*************************************************
[fv,fx] = Intial_Con(LevX,LevV,Deg,Lmax,Vmax,pde,...
    FMWT_COMP_x,FMWT_COMP_v);

%=============================================================
%% Step 2. Generate Sparse Grids/Hash Table
% Input:    Lev-Level of the Mesh
%           Deg-Degree of Polynomial
%           Dim-dimensionality
% Output: HASH and HashInv
%=============================================================
[HASH,HashInv] = HashTable(LevV,LevX,Deg,Dim);

% Generate the Initial Condition w.r.t Grids
fval = fv(HashInv.x1).*fx(HashInv.x2);

clear fv fx


%=============================================================
%% Step 3. Generate time-independent coefficient matrices
% Vlasolv Solver:
%   Operators:  vMassV: int_v v*l_i(v)*l_j(v)dv
%               GradV: int_v (l_i(v))'*l_j(v)dv
%               GradX: int_x (m_i(x))'*m_j(x)dx
% Poisson Solver:
%   Operators: DelaX: int_x (m_i(x))''*m_j(x)dx
% Maxwell Solver:
%   Operators: CurlCurlX: int_x curl(m_i(x))*m_j(x)dx
% Input: LevX, LevV, k, dim, Lmax, Vmax
% Output: 2D Matrices--vMassV,GradV,GradX,DeltaX, and CurlCurl(missing)
%=============================================================
[vMassV,GradV,GradX,DeltaX] = matrix_coeff_TI(LevX,LevV,Deg,Lmax,Vmax,...
    FMWT_COMP_x,FMWT_COMP_v);

% Generate A_encode for Time-independent Matrix
A_encode = GlobalMatrixSG(vMassV,GradX,HASH);

%====================================================================
%% Step 4. Generate time-independent global Matrix
% Compute the global matrix for spacial variables "x" by
%
% Poisson Solver: A_Poisson (Hash, Dim_x,k,LevX,DeltaX)
% or Maxwell Solver: A_Maxwell (Hash,Dim_x,k,LevX,nu,eps,CurlCurlX)
% Input: Hash, Dim_x,k,LevX,DeltaX,or nu, eps, CurlCurlX
% Output: A_Poisson or A_Maxwell
% Another Idea is to solve Poisson Equation on the finest full grid
%====================================================================
if DimX>1
    % Construct DeltaX for DimX
else
    A_Poisson = DeltaX;
end

%=============================================================
%% Step 5. Time Loop
%	Step 5.1 Vlasov Equation
%       Generate time dependent coefficient matrix
%       Generate global matrix A_Vlasov(Hash,coef_mat,Dim)
%       Apply A_Vlasov->f by RK
%	Step 5.2 Poisson Equation
%       Solve Poisson Equation sol_Poisson by A_Poisson(f)
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
[Meval_v,v_node,Meval_x,x_node]=matrix_plot(LevX,LevV,Deg,Lmax,Vmax,...
    FMWT_COMP_x,FMWT_COMP_v);
%---------------------
% plot for validating
%---------------------
[xx,vv]=meshgrid(x_node,v_node);
tmp=Multi_2D(Meval_v,Meval_x,fval,HASH,HashInv);

figure(1000)
mesh(xx,vv,reshape(tmp,Deg*2^LevX,Deg*2^LevV)','FaceColor','interp','EdgeColor','interp');
axis([0 Lmax -Vmax Vmax])
view(0,90)
colorbar

count=1;
% return

for L = 1:floor(TEND/dt)
    %=============================================================
    %% Step 5.3. Generate time-dependent coefficient matrix
    % Vlasolv Solver:
    %   Operators:  EMassX: int_x E(x)*m_i(x)*m_j(x)dx
    % Input: Lev, Deg, Lmax, Vmax
    % Output: EMassX
    % Note: E is solved by Poisson or Maxwell's equation
    %=============================================================
    % Poisson Solver: Solve E from 1-rho=1-int f dv
    [E,u] = PoissonSolve(LevX,Deg,Lmax,fval,A_Poisson,FMWT_COMP_x,Vmax);

    % Generate EMassX matrix
    EMassX = matrix_coeff_TD(LevX,Deg,Lmax,E,FMWT_COMP_x);
    
    % B_encode for Time-Dependent Matrices
    B_encode = GlobalMatrixSG(GradV,EMassX,HASH);
    C_encode=[A_encode B_encode]; % This step is GlobalVlasov
    
    %====================================
    % RK Time Stepping Method
    %====================================
    fval = TimeAdvance(C_encode,fval, dt);
 
    time(count) = L*dt;
    count=count+1;
    
   %---------------------
    % plot for validating
    %---------------------
    figure(1000)

    tmp=Multi_2D(Meval_v,Meval_x,fval,HASH,HashInv);
    mesh(xx,vv,reshape(tmp,Deg*2^LevX,Deg*2^LevV)','FaceColor','interp','EdgeColor','interp');
    axis([0 Lmax -Vmax Vmax])
    view(0,90)
    colorbar
    
    title(['Time at ',num2str(L*dt)])
    pause (0.01)
   
end
