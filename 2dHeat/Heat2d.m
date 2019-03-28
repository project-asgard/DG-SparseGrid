function pde = Heat2d
% Example PDE using the 2D (1x-1v) Heat Equation. This example PDE is
% time dependent (although not all the terms are time dependent). This
% implies the need for an initial condition. 

% df/dt = Delta f = d^2f/dx^2 + d^2f/dy^2
% First we set the dimensionality as x and v
% 
%% Setup the dimensions
% 
% Here we setup a 2D problem (x,v)

dim_v.name = 'v';
dim_v.BCL = 0;
dim_v.BCR = 0;
dim_v.domainMin = -1;
dim_v.domainMax = +1;
dim_v.lev = 2;
dim_v.deg = 2;
dim_v.FMWT = []; % Gets filled in later
dim_v.init_cond_fn = @Fv_0;

dim_x.name = 'x';
dim_x.BCL = 0;
dim_x.BCR = 0;
dim_x.domainMin = -1;
dim_x.domainMax = +1;
dim_x.lev = 2;
dim_x.deg = 2;
dim_x.FMWT = []; % Gets filled in later
dim_x.init_cond_fn = @Fx_0;

%%
% Add dimensions to the pde object
% Note that the order of the dimensions must be consistent with this across
% the remainder of this PDE.

pde.dimensions = {dim_v, dim_x}; % Order chosen here to match the old hard wired version

%% Setup the terms of the PDE
%
% Here we have 2 terms, with each term having nDims (x and v) operators.

%% 
% Setup the v.d_dx (v.MassV . GradX) term

term2_x.type = 1; % grad (see coeff_matrix.m for available types)
term2_x.G = @(x,t,dat) x*0+1; % G function for use in coeff_matrix construction.
term2_x.TD = 0; % Time dependent term or not.
term2_x.dat = []; % These are to be filled within the workflow for now
term2_x.LF = 0; % Use Lax-Friedrichs flux or not TODO : what should this value be?
term2_x.name = 'd_dx';

term2_v.type = 2; % mass (see coeff_matrix.m for available types)
term2_v.G = @(v,t,dat) v; % G function for use in coeff_matrix construction.
term2_v.TD = 0; % Time dependent term or not.
term2_v.dat = []; % These are to be filled within the workflow for now
term2_v.LF = 0; % Use Lax-Friedrichs flux or not.
term2_v.name = 'v';

term2 = {term2_v, term2_x};

%% 
% Setup the E.d_dv (E.MassX . GradV) term

term3_x.type = 2;
term3_x.G = @(x,t,dat) dat;
term3_x.TD = 1;
term3_x.dat = [];
term3_x.LF = 0;
term3_x.name = 'E';

term3_v.type = 1;
term3_v.G = @(v,t,dat) v*0+1;
term3_v.TD = 0;
term3_v.dat = [];
term3_v.LF = 0;
term3_v.name = 'd_dv';

term3 = {term3_v, term3_x};

%%
% Add terms to the pde object

pde.terms = {term2, term3};

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
analytic_v = @ExactFv;
analytic_t = @ExactFt;

pde.analytic_solutions_1D = {analytic_v,analytic_x,analytic_t};
pde.analytic_solution = @ExactF;

%% Other workflow options that should perhpas not be in the PDE?

pde.set_dt = @set_dt; % Function which accepts the pde (after being updated with CMD args).
pde.Ex = @Ex; % These can actually get absorbed into the G functions above.
pde.Et = @Et; % but I've not done it yet. 
pde.solvePoisson = 0; % Controls the "workflow" ... something we still don't know how to do generally. 
pde.applySpecifiedE = 1; % Controls the "workflow" ... something we still don't know how to do generally. 
pde.implicit = 0; % Can likely be removed and be a runtime argument. 
pde.checkAnalytic = 1; % Will only work if an analytic solution is provided within the PDE.

end

%% Define the various input functions specified above. 

function f=Fx_0(x,p)
% Initial condition for x variable
f=x.*0;
end
function f=Fv_0(v,p)
% Initial condition for v variable
f=v.*0;
end
function f=Fxv_0(x,v,p)
f=Fv_0(v).*Fx_0(x);
end

% Apply this specific E field
function f=Ex(x, p)
f=cos(pi*x);
end
function f=Et(t,p)
f=cos(t);
end
function f=E(x,t,p)
f=Ex(x).*Et(t);
end

%%
% Source terms are composed of fully seperable functions
% Source = source1 + source2 + source3

%%
% Source term 1
function f = source1t(t)
f = cos(t);
end
function f = source1x(x,p)
f = sin(pi*x);
end
function f = source1v(v,p)
f = sin(pi*v/5);
end
function f = source1(x,v,t)
f = source1x(x).*source1v(v).*source1t(t);
end

%%
% Source term 2
function f = source2t(t)
f = sin(t);
end
function f = source2x(x,p)
f = cos(pi*x);
end
function f = source2v(v,p)
f = pi*v.*sin(pi*v/5);
end
function f = source2(x,v,t)
f = source2x(x).*source2v(v).*source2t(t);
end

%%
% Source term 3
function f = source3t(t)
f = cos(t).*sin(t);
end
function f = source3x(x,p)
f = cos(pi*x).*sin(pi*x);
end
function f = source3v(v,p)
f = 1/5*pi*cos(pi*v/5);
end
function f = source3(x,v,t)
f = source3x(x).*source3v(v).*source3t(t);
end

%%
% Analytic Solution functions

function f=ExactFt(t)
f=sin(t);
end
function f=ExactFx(x,p)
f = sin(pi*x);
end
function f=ExactFv(v,p)
f = sin(pi*v/5);
end
function f=ExactF(x,v,t)
f = ExactFx(x).*ExactFv(v).*ExactFt(t);
end

%%
% Function to set time step
function dt=set_dt(pde)

Vmax = pde.dimensions{1}.domainMax;
Lmax = pde.dimensions{2}.domainMax;
LevX = pde.dimensions{2}.lev;
CFL = pde.CFL;
Deg = pde.dimensions{1}.deg;

dt = Lmax/2^LevX/Vmax/(2*Deg+1)*CFL;
end
