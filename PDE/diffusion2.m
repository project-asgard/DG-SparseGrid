function pde = diffusion2
% Example PDE using the 2D (1x-1y) Heat Equation. This example PDE is
% time dependent (although not all the terms are time dependent). This
% implies the need for an initial condition. 
% PDE:
% 
% df/dt = d^2 f/dx^2 + d^2 f/dy^2
%
% Domain is [0,1]x[0,1]
% Dirichlet boundary condition 
%
% Diffusion terms are dealt with via LDG, i.e., splitting into two first
% order equations:
%
% d^2 f / dx^2 becomes
%
% dq/dx with free (homogeneous Neumann BC)
%
% and
%
% q=df/dx with Dirichlet BCs specified by analytic solution
%
% Run with
%
% explicit
% fk6d(diffusion2,3,2,0.001);
%
% implicit
% fk6d(diffusion2,4,2,0.05,[],[],1,[],[],1.9);

pde.CFL = 0.01;

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
dim_x.init_cond_fn = @(x,p) soln_x(x)*soln_t(0);

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
dim_y.init_cond_fn = @(y,p) soln_y(y)*soln_t(0);

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

term1_x.type = 'diff';
% eq1 : dq/dx (flux equation)
term1_x.G1 = @(x,p,t,dat) x*0+1;
term1_x.LF1 = +1; % upwind right
term1_x.BCL1 = 'N';
term1_x.BCR1 = 'N';
% eq2 : df/dx (actual variable eqn)
term1_x.G2 = @(x,p,t,dat) x*0+1;
term1_x.LF2 = -1; % upwind left
term1_x.BCL2 = 'D';
term1_x.BCR2 = 'D';
term1_x.BCL2_fList = BCL_fList;
term1_x.BCR2_fList = BCR_fList;

term1 = term_fill({term1_x,[]});

%% 
% Setup the d^2_dy^2 term

term2_y.type = 'diff';
% eq1 : dq/dy (flux equation)
term2_y.G1 = @(y,p,t,dat) y*0+1;
term2_y.LF1 = +1; % upwind right
term2_y.BCL1 = 'N';
term2_y.BCR1 = 'N';
% eq2 : df/dy (actual variable eqn)
term2_y.G2 = @(y,p,t,dat) y*0+1;
term2_y.LF2 = -1; % upwind left
term2_y.BCL2 = 'D';
term2_y.BCR2 = 'D';
term2_y.BCL2_fList = BCL_fList;
term2_y.BCR2_fList = BCR_fList;

term2 = term_fill({[],term2_y});

%%
% Add terms to the pde object

 pde.terms = {term1, term2};

%% Construct some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

pde.params = params;

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
