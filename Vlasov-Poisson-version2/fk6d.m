function [fval,err] = fk6d(pde,Lev,Deg,TEND,quiet,compression)

%% MATLAB (reference) version of the DG-SG solver
% The default execution solves the Vlasov-Poisson system of equations
%
% $$f_t + v\frac{\partial f}{\partial x} + E\left(x,t\right)\frac{\partial
% f}{\partial v} f=0$$
%
% $$-\frac{\partial^2\phi}{\partial x^2} = \rho - 1$$
%
% $$E\left(x,t\right)= \frac{\partial\phi}{\partial x}$$
%
% where $\rho\left(x,t\right)=\int_v f(x,v,t) dv$.

format short e
addpath(genpath(pwd))

%% Step 1. Set input parameters
% pde :  A structure containing the initial condition and domain
% information. See PDE/vlasov4.m for and example. Note that this does not
% actually describe the PDE (that's done by the coefficient matricies), so
% we should probably rename this.
%
% Lev: Maximum level of the mesh in all dimensions.
%
% Deg: Degree of basis functions.
%
% TEND : End time of the simulation, i.e., run from t=0 to t=TEND in steps
% of dt.
%
% quiet : Print debugging statements or not.
%
% Compression : Choice of approach to constructing the A system matrix.

if ~exist('pde','var') || isempty(pde)
    % Equation setup
    pde = Vlasov4;
end
if ~exist('TEND','var') || isempty(TEND)
    % End time
    TEND = 1;
end
if ~exist('Lev','var') || isempty(Lev)
    % Number of levels
    Lev = 3;
end
if ~exist('Deg','var') || isempty(Deg)
    % Polynomial degree
    Deg = 2; % Deg = 2 Means Linear Element
end
if ~exist('quiet','var') || isempty(quiet)
    % Enable / disable print statements
    quiet = 0;
end
if ~exist('compression','var') || isempty(compression)
    % Use or not the compression reference version
    compression = 4;
end

% Get x and v domain ranges.
Vmax = pde.params.Vmax;
Lmax = pde.params.Lmax;

% Level information.
LevX = Lev;
LevV = Lev;
pde.params.LevX = LevX;
pde.params.LevV = LevV;

% Dimensionality.
Dim = 2;
DimX = 1;
DimV = 1;
pde.params.DimX = DimX;
pde.params.DimV = DimV;

% Degree
pde.params.Deg = Deg;

% Time step.
dt = Lmax/2^LevX/Vmax/(2*Deg+1);

%% Step 1.1. Setup the multi-wavelet transform in 1D (for each dimension).
%
% Input:  Deg and Lev
%
% Output: Convert Matrix FMWT_COMP

FMWT_COMP_x = OperatorTwoScale(Deg,2^LevX);
FMWT_COMP_v = OperatorTwoScale(Deg,2^LevV);


%% Step 1.2. Apply the mulit-wavelet transform to the initial conditions in each dimension.
% Generate the 1D initial conditions. Input: LevX,LevV,Deg,Lmax,Vmax,pde
% Output: fval (fv and fx)--intial condition f(x,v,t=0)
%         rho--intial condition rho(x,t=0)

if ~quiet; disp('[1.2] Setting up 1D initial conditions'); end
[fv,fx] = Intial_Con(LevX,LevV,Deg,Lmax,Vmax,pde,FMWT_COMP_x,FMWT_COMP_v);


%% Step 2. Generate Sparse-Grid (as the Hash + Connectivity tables).

%%% Construct forward and inverse hash tables.
if ~quiet; disp('[2.1] Constructing hash and inverse hash tables'); end
[HASH,HASHInv] = HashTable(Lev,Dim);
nHash = numel(HASHInv);

%%% Construct the connectivity.
if ~quiet; disp('[2.2] Constructing connectivity table'); end
Con2D = Connect2D(Lev,HASH,HASHInv);

%%% Get the multi-wavelet coefficient representation on the sparse-grid,
%%% i.e., above we transformed each of the initial condition 1D
%%% dependencies into the multi-wavelet basis. Here we combine those N 1D
%%% basis representations into the multi-D (here N=2, so 2D) initial
%%% condition, multi-wavelet representation of the initial condition
%%% specified in PDE.
if ~quiet; disp('[2.3] Calculate 2D initial condition on the sparse-grid'); end
fval = initial_condition_vector(fx,fv,Deg,Dim,HASHInv);

clear fv fx

%% Step 3. Generate time-independent coefficient matrices
% Vlasolv Solver:
%   Operators:
%               vMassV: int_v v*l_i(v)*l_j(v)dv GradV: int_v
%               (l_i(v))'*l_j(v)dv GradX: int_x (m_i(x))'*m_j(x)dx
% Poisson Solver:
%               Operators: DelaX: int_x (m_i(x))''*m_j(x)dx
% Input:
%               LevX, LevV, k, dim, Lmax, Vmax
% Output:
%               2D Matrices--vMassV,GradV,GradX,DeltaX

%%% Build the time independent coefficient matricies.
if ~quiet; disp('[3.1] Calculate time independent matrix coefficients'); end
[vMassV,GradV,GradX,DeltaX] = matrix_coeff_TI(LevX,LevV,Deg,Lmax,Vmax,...
    FMWT_COMP_x,FMWT_COMP_v);

%%% Generate A_encode / A_data time independent data structures.
if ~quiet; disp('[3.2] Generate A_encode data structure for time independent coefficients'); end
if compression == 3
    A_encode=GlobalMatrixSG(vMassV,GradX,HASHInv,Con2D,Deg);
else
    % A_data is constructed only once per grid refinement, so can be done
    % on the host side.
    A_data = GlobalMatrixSG_SlowVersion(HASHInv,Con2D,Deg,compression);
end

%% Step 4. Generate time-independent global matrix
% Compute the global matrix for spatial variables "x" by
%
% Poisson Solver: A_Poisson (Hash, Dim_x,k,LevX,DeltaX)
%
% Input: Hash, Dim_x,k,LevX,DeltaX,or nu, eps, CurlCurlX Output: A_Poisson
% Another Idea is to solve Poisson Equation on the finest full grid

if ~quiet; disp('[4] Construct matrix for Poisson solve'); end
if DimX>1
    % Construct DeltaX for DimX
else
    A_Poisson = DeltaX;
end


%% Step 5. Time Loop
%	Step 5.1 Vlasov Equation
%       Generate time dependent coefficient matrix Generate global matrix
%       A_Vlasov(Hash,coef_mat,Dim) Apply A_Vlasov->f by RK
%	Step 5.2 Poisson Equation
%       Solve Poisson Equation sol_Poisson by A_Poisson(f) Compute
%       E=(sol_Poisson)'


% At time = 0 plotting.
if ~quiet; disp('[5.0] Plotting intial condition'); end

if ~quiet
    % Prepare plotting data.
    [Meval_v,v_node,Meval_x,x_node]=matrix_plot(LevX,LevV,Deg,Lmax,Vmax,...
        FMWT_COMP_x,FMWT_COMP_v);
    
    % Plot initial condition.
    [xx,vv]=meshgrid(x_node,v_node);
    tmp=Multi_2D(Meval_v,Meval_x,fval,HASHInv,Lev,Deg);
    
    figure(1000)
    mesh(xx,vv,reshape(tmp,Deg*2^LevX,Deg*2^LevV)','FaceColor','interp','EdgeColor','interp');
    axis([0 Lmax -Vmax Vmax])
    view(0,90)
    colorbar
end


% Write the initial condition to file.
write_fval = 0;
if write_fval; write_fval_to_file(fval,Lev,Deg,0); end

count=1;
plotFreq = 10;

if ~quiet; disp('[7] Advancing time ...'); end
for L = 1:floor(TEND/dt)
    
    time(count) = L*dt;
    timeStr = sprintf('Step %i of %i',L,floor(TEND/dt));
    
    if ~quiet; disp(timeStr); end
    
    if pde.solvePoisson
        %%% Solve Poisson to get E (from 1-rho=1-int f dv)
        if ~quiet; disp('    [a] Solve poisson to get E'); end
        E = PoissonSolve(LevX,Deg,Lmax,fval,A_Poisson,FMWT_COMP_x,Vmax);
    end
    
    if pde.applySpecifiedE
        %%% Apply specified E
        if ~quiet; disp('    [a] Apply specified E'); end
        % Does this take E(x) -> E(lev,pos,deg) ?
        E = ProjCoef2Wav_v2(LevX,Deg,0,Lmax,pde.exactE);
    end
    
    %%% Generate EMassX time dependent coefficient matrix.
    if ~quiet; disp('    [b] Calculate time dependent matrix coeffs'); end
    EMassX = matrix_coeff_TD(LevX,Deg,Lmax,E,FMWT_COMP_x);
    
    %%% Update A_encode for time-dependent coefficient matricies.
    if ~quiet; disp('    [c] Generate A_encode for time-dependent coeffs'); end
    if compression == 3
        B_encode = GlobalMatrixSG(GradV,EMassX,HASHInv,Con2D,Deg);
        C_encode=[A_encode B_encode];
    else
        
    end
    
    %%% Advance Vlasov in time with RK3 time stepping method.
    if ~quiet; disp('    [d] RK3 time step'); end
    if compression == 3
        fval = TimeAdvance(C_encode,fval,time(count),dt,compression,Deg,pde,HASHInv);
    else
        A_data.vMassV    = vMassV;
        A_data.GradX     = GradX;
        A_data.GradV     = GradV;
        A_data.EMassX    = EMassX;
        
        % Write the A_data structure components for use in HPC version.
        write_A_data = 0;
        if write_A_data && L==1; write_A_data_to_file(A_data,Lev,Deg); end
        
        fval = TimeAdvance(A_data,fval,time(count),dt,compression,Deg,pde,HASHInv);
        
    end
    
    % Write the present fval to file.
    if write_fval; write_fval_to_file(fval,Lev,Deg,L); end
    
    
    %%% Plot results
    if mod(L,plotFreq)==0 && ~quiet
        
        figure(1000)
        
        tmp=Multi_2D(Meval_v,Meval_x,fval,HASHInv,Lev,Deg);
        mesh(xx,vv,reshape(tmp,Deg*2^LevX,Deg*2^LevV)','FaceColor','interp','EdgeColor','interp');
        axis([0 Lmax -Vmax Vmax])
        view(0,90)
        colorbar
        
        title(['Time at ', timeStr])
        pause (0.01)
    end
    
    if pde.checkAnalytic
        %%% Check the solution with the analytic solution
        fval_analytic = source_vector2(LevX,LevV,Deg,HASHInv,pde,time(count));
        err = sqrt(mean(((fval - fval_analytic)/fval_analytic).^2));
    end
    
    count=count+1;
  
end

end

