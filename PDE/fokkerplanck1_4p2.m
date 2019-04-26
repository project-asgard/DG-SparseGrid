function pde = fokkerplanck1_4p2
% 1D pitch angle collisional term 
% df/dt == d/dz ( (1-z^2) df/dz ) 
% 
% Here we use LDG for this second order system. We impose homogeneous
% Neumann BCs on f
%
% d/dz( (1-z^2) df/dz ) becomes
%
% d/dz (1-z^2)*q  with free (homogeneous Neumann BC)
%
% and the flux is 
%
% q=df/fz  with homogeneous Dirichlet BC
%
% Run with
%
% explicit
% fk6d(fokkerplanck1_4p2,4,4,0.01);
%
% implicit
% fk6d(fokkerplanck1_4p2,5,4,0.01,[],[],1,[],[],0.5);

pde.CFL = 0.005;


%% Setup the dimensions
% 
% Here we setup a 1D problem (x)

    function ret = soln(z,t)
        
        h = [3,0.5,1,0.7,3,0,3];
        
        ret = zeros(size(z));
        for l=1:numel(h)
            
            L = l-1;
            P_m = legendre(L,z); % Use matlab rather than Lin's legendre.
            P = P_m(1,:)';
            
            ret = ret + h(l) * P * exp(-L*(L+1)*t);
            
        end
        
    end

dim_z.BCL = 'D'; % dirichlet
dim_z.BCR = 'D';
dim_z.domainMin = -1;
dim_z.domainMax = +1;
dim_z.init_cond_fn = @(z,p) soln(z,0);

%%
% Add dimensions to the pde object
% Note that the order of the dimensions must be consistent with this across
% the remainder of this PDE.

pde.dimensions = {dim_z};

%% Setup the terms of the PDE
%
% Here we have 1 term1, having only nDims=1 (x) operators.

%% 
% d/dz( (1-z^2) df/dz )

termC_z.type = 'diff';
% eq1 : 1 * d/dx (1-z^2) q
termC_z.G1 = @(z,p,t,dat) 1-z.^2;
termC_z.LF1 = -1; % upwind left
termC_z.BCL1 = 'D';
termC_z.BCR1 = 'D';
% eq2 : q = df/dx 
termC_z.G2 = @(z,p,t,dat) z*0+1;
termC_z.LF2 = +1; % upwind right
termC_z.BCL2 = 'N';
termC_z.BCR2 = 'N';

termC = term_fill({termC_z});

%%
% Add terms to the pde object

pde.terms = {termC};

%% Construct some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

pde.params = params;

%% Add an arbitrary number of sources to the RHS of the PDE
% Each source term must have nDims + 1

%%
% Add sources to the pde data structure
pde.sources = {};

%% Define the analytic solution (optional).
% This requires nDims+time function handles.

pde.analytic_solutions_1D = { ...
    @(z,p,t) soln(z,t), ...
    @(t,p) 1 
    };

%%
% Function to set time step

    function dt=set_dt(pde)
        dims = pde.dimensions;
        % for Diffusion equation: dt = C * dx^2
        lev = dims{1}.lev;
        xMax = dims{1}.domainMax;
        xMin = dims{1}.domainMin;
        xRange = xMax-xMin;
        CFL = pde.CFL;
        dx = xRange/2^lev;
        dt = CFL*dx^2;
    end

pde.set_dt = @set_dt;

end

