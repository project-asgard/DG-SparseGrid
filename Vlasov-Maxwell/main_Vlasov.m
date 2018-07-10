%**************************************************************************
% This MATLAB-Interface for DG-SG Vlasov Equations
% is an example to solve Vlasov Equations with Poisson/Maxwell Equations
% Vlasov Equation with Poisson solver:
% (1). f_t + v*grad_x f + E*grad_v f=0;
% (2). -Del (phi) = 1 - int_v f dv and E= (phi)'
% Vlasov Equation with Maxwell solver:
% (1). f_t + v*grad_x f + (E+B)*grad_v f=0;
% (2). Maxwell Eq.
% Input: 
%           Solver = 'VM' -- Vlasov Maxwell Solver
%                   'VP' -- Vlasov Poisson Solver
%           pde = choosing your test from folder "PDE"
%           DimX = (Dim_x1,Dim_x2,Dim_x3)
%           DimV = (Dim_v1,Dim_v2,Dim_v3)
% The documentation is located at ./Doc/html/Introduction.html
%**************************************************************************
clc
clear
close all
format short e

addpath(genpath(pwd))

%------- Input Parameters
Solver = 'VP';
pde = VlasovPoisson4;

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
%**************************
% Test PDE and Ending Time
%**************************
if Solver == 'VM' 
    if isempty('pde')==0
        pde = VlasovMaxwell1;
    end
    Vmax = pde.Vmax;
    Lmax = pde.Lmax;
    TEND = 10;
    
elseif Solver == 'VP' 
    if isempty('pde')==0
        pde = VlasovPoisson4;
    end
    Vmax = pde.Vmax;
    Lmax = pde.Lmax;
    TEND = 10;
    Cb1=0;Cb2=0;Cb3=0;

end

%*********************
% For dimensionality
%*********************
if pde.DimV == 1
    Cv2=0;Cv3=0;
end
if (pde.DimX ==1) 
    if (Solver == 'VP')
    Cx2=0;Cx3=0;

    elseif (Solver == 'VM')
    Cx1=0;Cx3=0;
    end
end

% Dimensionality
Dim = pde.DimX+pde.DimV;

%*********************
% Level information
%*********************
Lev = 3;

%*********************
% Polynomial Degree
%*********************
Deg = 2;% Deg = 2 Means Linear Element

%*************************************
% Time step
% Need future effort for CFL condition
%*************************************
dt = Lmax/2^Lev/Vmax/(2*Deg+1);

%*************************************************
%% Step 1.1. Set Up Matrices for Multi-wavelet
% Input:  Deg and Lev
% Output: Convert Matrix FMWT_COMP
%*************************************************
FMWT_COMP_x = OperatorTwoScale(Deg,2^Lev);
FMWT_COMP_v = OperatorTwoScale(Deg,2^Lev);

%*************************************************
%% Step 1.2. Initial Data according to Hash
% Generate the initial condition
% Input: LevX,LevV,Deg,Lmax,Vmax,pde
% Output: fval (fv and fx)--intial condition f(x,v,t=0)
%         rho--intial condition rho(x,t=0)
%*************************************************
InCond = Intial_Con(Lev,Deg,Lmax,Vmax,pde,...
    FMWT_COMP_x,FMWT_COMP_v,Solver);

%=============================================================
%% Step 2. Generate Sparse Grids/Hash Table
% Input:    Lev-Level of the Mesh
%           Deg-Degree of Polynomial
%           Dim-dimensionality
% Output: HASH and HashInv
%=============================================================
[HASH,inverseHash] = HashTable(Lev,Dim);

% Generate the Initial Condition w.r.t Grids
% needs some work over here!!
%%fval = fv(HashInv.x1).*fx(HashInv.x2);
% % finit = Initial_Con_Grid(HASH,inverseHash,InCond);
% % clear fv fx


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
% Output: 2D Matrices--vMassV,GradV,GradX,DeltaX
%=============================================================
% [vMassV,GradV,GradX,DeltaX] = matrix_coeff_TI(Lev,Deg,Lmax,Vmax,...
%     FMWT_COMP_x,FMWT_COMP_v);
Matrix_TI=matrix_coeff_TI_v2(Lev,Deg,Lmax,Vmax,FMWT_COMP_x,FMWT_COMP_v);

% for generating 6D matrix
A_encode = GlobalMatrixSG(I1,I2,I3,I4,I5,I6,Hash);
return
% Generate A_encode for Time-independent Matrix
A_encode = GlobalMatrixSG(vMassV,GradX,HASH);

% Building 12 matrices corresponding to each terms indicated as notes
% Note: df/dxj and df/dvj mean the partial derivatives
% I2 =  (v1*(df/dx1),w)
% I3 =  (v2*(df/dx2),w)
% I4 =  (v3*(df/dx3),w)
% I5 =  (E1*(df/dv1),w)
% I6 =  (v2*B3*(df/dv1),w)
% I7 = -(v3*B2*(df/dv1),w)
% I8 =  (E2*(df/dv2),w)
% I9 =  (v3*B1*(df/dv2),w)
% I10 = -(v1*B3*(df/dv2),w)
% I11 = (E3*(df/dv3),w)
% I12 = (v1*B2*(df/dv3),w)
% I13 = -(v2*B1*(df/dv3),w)

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
if Solver == 'VM'
    A_Poisson = MaxwellSolver(DeltaX,DimX);
elseif Solver =='VP'
    A_Poisson = PoissonSolver(DeltaX,DimX);
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
    E = PoissonSolve(LevX,Deg,Lmax,fval,A_Poisson,FMWT_COMP_x,Vmax);

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
