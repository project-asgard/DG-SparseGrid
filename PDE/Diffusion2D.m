function pde = Diffusion2D
% Example PDE using the 2D (1x-1y) Heat Equation. This example PDE is
% time dependent (although not all the terms are time dependent). This
% implies the need for an initial condition. 
% PDE: df/dt = d^2 f/dx^2 + d^2 f/dy^2
% Domain is [0,1]x[0,1]
% Dirichlet boundary condition

%% Setup the dimensions
% 
% Here we setup a 2D problem (x,v)
lev = 2;
deg = 2;

dim_x.name = 'x';
dim_x.BCL = 1; % Dirichlet
dim_x.BCR = 1; % Dirichlet
dim_x.domainMin = 0;
dim_x.domainMax = 1;
dim_x.lev = lev;
dim_x.deg = lev;
dim_x.FMWT = []; % Gets filled in later
dim_x.init_cond_fn = @Fx_0;

dim_y.name = 'y';
dim_y.BCL = 1;
dim_y.BCR = 1;
dim_y.domainMin = 0;
dim_y.domainMax = 1;
dim_y.lev = lev;
dim_y.deg = deg;
dim_y.FMWT = []; % Gets filled in later
dim_y.init_cond_fn = @Fy_0;

%%
% Add dimensions to the pde object
% Note that the order of the dimensions must be consistent with this across
% the remainder of this PDE.

pde.dimensions = {dim_x, dim_y}; % Order chosen here to match the old hard wired version

%% Setup the terms of the PDE
%
% Here we have 1 term, with each term having nDims (x and y) operators.

%% 
% Setup the v.d_dx (v.MassV . GradX) term
term_1D.dat = [];
term_1D.LF = 1;       % Upwind Flux
term_1D.G = @(x,t,y)1; % Delta Operator 
term_1D.type = 3;      % Delta Operator ::  Let this denote the derivative order
pde.term_1D = term_1D;

% term2_x.type = 1; % grad (see coeff_matrix.m for available types)
% term2_x.G = @(x,t,dat) x*0+1; % G function for use in coeff_matrix construction.
% term2_x.TD = 0; % Time dependent term or not.
% term2_x.dat = []; % These are to be filled within the workflow for now
% term2_x.LF = 1; % Use Lax-Friedrichs flux or not TODO : what should this value be?
% term2_x.name = 'd_dx';
% 
% term2_y.type = 2; % mass (see coeff_matrix.m for available types)
% term2_y.G = @(v,t,dat) 1; % G function for use in coeff_matrix construction.
% term2_y.TD = 0; % Time dependent term or not.
% term2_y.dat = []; % These are to be filled within the workflow for now
% term2_y.LF = 0; % Use Lax-Friedrichs flux or not.
% term2_y.name = 'v';
% 
% term2 = {term2_y, term2_x};
% 
% %% 
% % Setup the E.d_dv (E.MassX . GradV) term
% 
% term3_x.type = 2;
% term3_x.G = @(x,t,dat) dat;
% term3_x.TD = 1;
% term3_x.dat = [];
% term3_x.LF = 0;
% term3_x.name = 'E';
% 
% term3_v.type = 1;
% term3_v.G = @(v,t,dat) v*0+1;
% term3_v.TD = 0;
% term3_v.dat = [];
% term3_v.LF = 0;
% term3_v.name = 'd_dv';
% 
% term3 = {term3_v, term3_x};

%%
% Add terms to the pde object

% pde.terms = {term2, term3};

%% Construct some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

pde.params = params;

%% Add an arbitrary number of sources to the RHS of the PDE
% Each source term must have nDims + 1 (which here is 2+1 / (x,v) + time) functions describing the
% variation of each source term with each dimension and time.
% Here we define 3 source terms.

%%
% Source 1
s1x = @source1x;
s1v = @source1v;
s1t = @source1t;
source1 = {s1v,s1x,s1t};

%%
% Source 2
s2x = @source2x;
s2v = @source2v;
s2t = @source2t;
source2 = {s2v,s2x,s2t};

%%
% Source 3
s3x = @source3x;
s3v = @source3v;
s3t = @source3t;
source3 = {s3v,s3x,s3t};

%%
% Add sources to the pde data structure
pde.sources = {source1,source2,source3};

%% Define the analytic solution (optional).
% This requires nDims+time function handles.

analytic_x = @ExactFx;
analytic_y = @ExactFy;
analytic_t = @ExactFt;

pde.analytic_solutions_1D = {analytic_y,analytic_x,analytic_t};
pde.analytic_solution = @ExactF;

%% Other workflow options that should perhpas not be in the PDE?
% Need some work here
pde.set_dt = @set_dt; % Function which accepts the pde (after being updated with CMD args).
pde.Ex = @Ex; % These can actually get absorbed into the G functions above.
pde.Et = @Et; % but I've not done it yet. 
pde.solvePoisson = 0; % Controls the "workflow" ... something we still don't know how to do generally. 
pde.applySpecifiedE = 1; % Controls the "workflow" ... something we still don't know how to do generally. 
pde.implicit = 0; % Can likely be removed and be a runtime argument. 
pde.checkAnalytic = 1; % Will only work if an analytic solution is provided within the PDE.

end

%% Define the various input functions specified above. 

function f = Fx_0(x,p)
% Initial condition for x variable
f = cos(pi*x);
end
function f = Fy_0(y,p)
% Initial condition for v variable
f = cos(pi*y);
end
function f = Fxy_0(x,y,p)
f = Fy_0(y).*Fx_0(x);
end


%%
% Analytic Solution functions
% f(x,y,t) = f(x)f(y)f(t)
function f=ExactFt(t,p)
f = exp(-2*pi^2*t);
end
function f=ExactFx(x,p)
f = cos(pi*x);
end
function f=ExactFy(y,p)
f = cos(pi*y);
end
function f=ExactF(x,y,t)
f = ExactFx(x).*ExactFy(y).*ExactFt(t);
end

%% Source term
% df/dt - d^2 f/dx^2 = 0 for this test
% But not true for other Manu-Sol
function f = sourcet(t)
f = t-t;
end
function f = sourcex(x)
f = x-x;
end
function f = sourcey(y)
f = y-y;
end
function f = source(x,v,t)
f = sourcex(x).*sourcey(v).*sourcet(t);
end
%%
% Function to set time step
function dt=set_dt(pde)

LXmax = pde.dimensions{1}.domainMax;
LYmax = pde.dimensions{2}.domainMax;
LevX = pde.dimensions{2}.lev;
CFL = pde.CFL;
Deg = pde.dimensions{1}.deg;

% for Diffusion equation: dt = C * dx^2
dt = LXmax/2^LevX/Vmax/(2*Deg+1)*CFL;
end
