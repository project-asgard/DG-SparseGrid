function pde = continuity2
% 2D test case using continuity equation, i.e.,
% df/dt + v_x * df/dx + v_y * df/dy == 0
%
% Run with ...
%  fk6d(continuity2,4,3,0.01,[],[],0,[]);

%% Setup the dimensions
%
% Here we setup a 2D problem (x,y)
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

dim_x.name = 'x';
dim_x.BCL = 'P'; % periodic
dim_x.BCR = 'P';
dim_x.domainMin = -1;
dim_x.domainMax = +1;
dim_x.lev = 2;
dim_x.deg = 2;
dim_x.FMWT = []; % Gets filled in later
dim_x.init_cond_fn = @(x,p) x.*0;

dim_x = checkDimension(dim_x);

dim_y.name = 'y';
dim_y.BCL = 'P'; % periodic
dim_y.BCR = 'P';
dim_y.domainMin = -2;
dim_y.domainMax = +2;
dim_y.lev = 2;
dim_y.deg = 2;
dim_y.FMWT = []; % Gets filled in later
dim_y.init_cond_fn = @(y,p) y.*0;

dim_y = checkDimension(dim_y);

%%
% Add dimensions to the pde object
% Note that the order of the dimensions must be consistent with this across
% the remainder of this PDE.

pde.dimensions = {dim_x,dim_y};

%% Setup the terms of the PDE
%
% Here we have 2 terms, having only nDims=2 (x,y) operators.

%%
% Setup the v_x * d_dx (v_x . GradX . MassY ) term

term2_x.type = 'grad'; % grad (see coeff_matrix.m for available types)
term2_x.G = @(x,p,t,dat) x*0-1; % G function for use in coeff_matrix construction.
term2_x.TD = 0; % Time dependent term or not.
term2_x.dat = []; % These are to be filled within the workflow for now
term2_x.LF = 0; % Use Lax-Friedrichs flux or not TODO : what should this value be?
term2_x.name = 'v_x.d_dx';

term2 = term_fill({term2_x,[]});


%%
% Setup the v_y * d_dy (v_y . MassX . GradY) term

term3_y.type = 'grad'; % grad (see coeff_matrix.m for available types)
term3_y.G = @(y,p,t,dat) y*0-1; % G function for use in coeff_matrix construction.
term3_y.TD = 0; % Time dependent term or not.
term3_y.dat = []; % These are to be filled within the workflow for now
term3_y.LF = 0; % Use Lax-Friedrichs flux or not TODO : what should this value be?
term3_y.name = 'v_y.d_dy';

term3 = term_fill({[],term3_y});

%%
% Add terms to the pde object

pde.terms = {term2,term3};

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
% Source term 1
source1 = { ...
    @(x,p) cos(pi*x),   ...   % s1x
    @(y,p) sin(2*pi*y), ...   % s1y
    @(t)   2*cos(2*t)   ...   % s1t
    };

%%
% Source term 2
source2 = { ...
    @(x,p)  cos(pi*x),    ...   % s2x
    @(y,p)  cos(2*pi*y),  ...   % s2y
    @(t)    2*pi*sin(2*t) ...   % s2t
    };

%%
% Source term 3
source3 = { ...
    @(x,p)  sin(pi*x),   ...  % s3x
    @(y,p)  sin(2*pi*y), ...  % s3y
    @(t)    -pi*sin(2*t) ...  % s3t
    };

%%
% Add sources to the pde data structure
pde.sources = {source1,source2,source3};

%% Define the analytic solution (optional).
% This requires nDims+time function handles.

pde.analytic_solutions_1D = { ...
    @(x,p,t)  cos(pi*x),   ... % a_x
    @(y,p,t)  sin(2*pi*y), ... % a_y
    @(t)    sin(2*t)     ... % a_t
    };

%% Other workflow options that should perhpas not be in the PDE?

pde.set_dt = @set_dt; % Function which accepts the pde (after being updated with CMD args).
pde.solvePoisson = 0; % Controls the "workflow" ... something we still don't know how to do generally.
pde.applySpecifiedE = 0; % Controls the "workflow" ... something we still don't know how to do generally.
pde.implicit = 0; % Can likely be removed and be a runtime argument.
pde.checkAnalytic = 1; % Will only work if an analytic solution is provided within the PDE.

end

%% Define the various input functions specified above.

% %%
% % Initial conditions for each dimension
% 
% function f=f0_x(x,p); f=x.*0; end
% function f=f0_y(y,p); f=y.*0; end

%%
% Function to set time step

function dt=set_dt(pde)

Lmax = pde.dimensions{1}.domainMax;
LevX = pde.dimensions{1}.lev;
CFL = pde.CFL;

dt = Lmax/2^LevX*CFL;
end
