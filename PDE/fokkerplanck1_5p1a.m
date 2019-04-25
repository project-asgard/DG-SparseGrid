function pde = fokkerplanck1_5p1a
% Test the momentum dynamics for the RE problem for E = R = 0
%
% x^2 * df/dt == d/dx * x^2 ( psi(x)/x * df/dx + 2*psi(x)*f )
%             == d/dx*x^2*psi(x)/x*df/dx  +  d/dx*x^2*2*psi(x)*f
%
% (note the x^2 moved to the left side)
%
%
% Run with
%
% explicit
% fk6d(fokkerplanck1_5p1a,5,3,0.01)
%
% implicit
% fk6d(fokkerplanck1_5p1a,5,4,3,[],[],1,'SG',[],1.5)

pde.CFL = 0.01;

%% Setup the dimensions
% 
% Here we setup a 1D problem (x)

    function ret = psi(x,t)
        
        phi = erf(x);
        dphi_dx = 2/sqrt(pi) * exp(-x.^2);

        ret = 1/(2*x.^2) * (phi - x.*dphi_dx)
    end

    function ret = f0(x)
        a = 2;
        ret = 4.0/(sqrt(pi)*a^3) * exp(-x.^2/a^2);
    end

    function ret = soln(x,t)
        ret = 4/sqrt(pi) * exp(-x.^2);
    end

dim_z.name = 'z';
dim_z.BCL = 'N'; % neumann
dim_z.BCR = 'N'; % not set (equivalent to neumann)
dim_z.domainMin = -1;
dim_z.domainMax = +1;
dim_z.lev = 2;
dim_z.deg = 2;
dim_z.FMWT = []; % Gets filled in later
dim_z.init_cond_fn = @(z,p) f0(z);

dim_z = checkDimension(dim_z);

%%
% Add dimensions to the pde object
% Note that the order of the dimensions must be consistent with this across
% the remainder of this PDE.

pde.dimensions = {dim_z};

%% Setup the terms of the PDE

%% 
% d/dx*x^2*psi(x)/x*df/dx

term1_z.type = 'diff'; % grad (see coeff_matrix.m for available types)
term1_z.G = @(x,p,t,dat) x.^2.*psi(x)./x; % G function for use in coeff_matrix construction.
term1_z.LF = -1; % Upwind 

term2 = term_fill({term1_z});

%%
% Add terms to the pde object

pde.terms = {term2};

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

%% Other workflow options that should perhpas not be in the PDE?

pde.set_dt = @set_dt; % Function which accepts the pde (after being updated with CMD args).
pde.solvePoisson = 0; % Controls the "workflow" ... something we still don't know how to do generally. 
pde.applySpecifiedE = 0; % Controls the "workflow" ... something we still don't know how to do generally. 
pde.implicit = 0; % Can likely be removed and be a runtime argument. 
pde.checkAnalytic = 1; % Will only work if an analytic solution is provided within the PDE.

end

%% Define the various input functions specified above. 

%%
% Function to set time step
function dt=set_dt(pde)

dims = pde.dimensions;
xRange = dims{1}.domainMax-dims{1}.domainMin;
lev = dims{1}.lev;
CFL = pde.CFL;
dx = xRange/2^lev;
dt = CFL * dx;

end
