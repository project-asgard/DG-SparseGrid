function pde = diffusion2
% Example PDE using the 2D (1x-1y) Heat Equation. This example PDE is
% time dependent (although not all the terms are time dependent). This
% implies the need for an initial condition. 
% PDE: df/dt = d^2 f/dx^2 + d^2 f/dy^2
% Domain is [0,1]x[0,1]
% Dirichlet boundary condition 
% ToDo: need some effort for naming, boundary conditions, source terms
%
% Run with
%
% explicit
% fk6d(diffusion2,3,2,0.001);
%
% implicit
% fk6d(diffusion2,4,2,0.05,[],[],1,[],[],1.9);

pde.CFL = 0.01;

lev = 5;
deg = 2;

%% Setup the dimensions
% 
% Here we setup a 2D problem (x,y)

soln_x = @(x) cos(pi*x);
soln_y = @(y) cos(pi*y);
soln_t = @(t) exp(-2*pi^2*t);

BCFunc = @(x) soln_x(x);
BCFunc_t = @(t) soln_t(t);

% Domain is (a,b)x(c,d)

% The function is defined for the plane
% x = a and x = b
BCL_fList = { ...
    @(x,p,t) BCFunc(x), ... % replace x by a
    @(y,p,t) BCFunc(y), ...
    @(t,p) BCFunc_t(t)
    };

BCR_fList = { ...
    @(x,p,t) BCFunc(x), ... % replace x by b
    @(y,p,t) BCFunc(y), ...
    @(t,p) BCFunc_t(t)
    };

dim_x.name = 'x';
dim_x.BCL = 'D'; % Dirichlet
dim_x.BCL_fList = BCL_fList;
dim_x.BCR = 'D'; % Dirichlet
dim_x.BCR_fList = BCR_fList;
dim_x.domainMin = 0;
dim_x.domainMax = 1;
dim_x.lev = lev;
dim_x.deg = deg;
dim_x.FMWT = []; % Gets filled in later
dim_x.init_cond_fn = @(x,p) soln_x(x)*soln_t(0);

dim_x = checkDimension(dim_x);

% The function is defined for the plane
% y = c and y = d
BCL_fList = { ...
    @(x,p,t) BCFunc(x), ...
    @(y,p,t) BCFunc(y), ... % replace y by c
    @(t,p) BCFunc_t(t)
    };

BCR_fList = { ...
    @(x,p,t) BCFunc(x), ...
    @(y,p,t) BCFunc(y), ...  % replace y by d
    @(t,p) BCFunc_t(t)
    };

dim_y.name = 'y';
dim_y.BCL = 'D';
dim_y.BCL_fList = BCL_fList;
dim_y.BCR = 'D';
dim_y.BCR_fList = BCR_fList;
dim_y.domainMin = 0;
dim_y.domainMax = 1;
dim_y.lev = lev;
dim_y.deg = deg;
dim_y.FMWT = []; % Gets filled in later
dim_y.init_cond_fn = @(y,p) soln_y(y)*soln_t(0);

dim_y = checkDimension(dim_y);

%%
% Add dimensions to the pde object
% Note that the order of the dimensions must be consistent with this across
% the remainder of this PDE.

pde.dimensions = {dim_x, dim_y};

%% Setup the terms of the PDE
%
% Here we have 1 term, with each term having nDims (x and y) operators.

%% 
% Setup the d^2_dx^2 term

term1_x.dat = [];
term1_x.LF = 0;         % Upwind Flux
term1_x.G = @(x,p,t,dat) x.*0 + 1; % Delta Operator 
term1_x.type = 'diff';       % Delta Operator ::  Let this denote the derivative order
term1_x.TD = 0;

term1 = term_fill({term1_x,[]});

%% 
% Setup the d^2_dy^2 term

term2_y.dat = [];
term2_y.LF = 0;         % Upwind Flux
term2_y.G = @(x,p,t,dat) x.*0 + 1; % Delta Operator 
term2_y.type = 'diff';        % Delta Operator ::  Let this denote the derivative order
term2_y.TD = 0;

term2 = term_fill({[],term2_y});

%%
% Add terms to the pde object

 pde.terms = {term1, term2};

%% Construct some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

pde.params = params;

%% Add an arbitrary number of sources to the RHS of the PDE
% Each source term must have nDims + 1 (which here is 2+1 / (x,v) + time) functions describing the
% variation of each source term with each dimension and time.
% Here we define 3 source terms.

% %%
% % Source 1
% s1x = @source1x;
% s1v = @source1v;
% s1t = @source1t;
% source1 = {s1v,s1x,s1t};


%%
% Add sources to the pde data structure
pde.sources = {};

%% Define the analytic solution (optional).
% This requires nDims+time function handles.

pde.analytic_solutions_1D = { ...
    @(x,p,t) soln_x(x), ...
    @(y,p,t) soln_y(y), ... 
    @(t,p) soln_t(t) 
    };

%% Other workflow options that should perhpas not be in the PDE?
% Need some work here
pde.set_dt = @set_dt; % Function which accepts the pde (after being updated with CMD args).
pde.Ex = @Ex; % These can actually get absorbed into the G functions above.
pde.Et = @Et; % but I've not done it yet. 
pde.solvePoisson = 0; % Controls the "workflow" ... something we still don't know how to do generally. 
pde.applySpecifiedE = 0; % Controls the "workflow" ... something we still don't know how to do generally. 
pde.implicit = 0; % Can likely be removed and be a runtime argument. 
pde.checkAnalytic = 1; % Will only work if an analytic solution is provided within the PDE.

end

%%
% Function to set time step
function dt=set_dt(pde)

dims = pde.dimensions;

% for Diffusion equation: dt = C * dx^2

lev = dims{1}.lev;
CFL = pde.CFL;
dx = 1/2^lev;
dt = CFL*dx^2;

end
