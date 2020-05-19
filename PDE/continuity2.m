function pde = continuity2
% 2D test case using continuity equation, i.e.,
%
% df/dt == -df/dx - df/dy
%
% Run with
%
% explicit
% asgard(continuity2,'lev',3,'deg',3)
%
% implicit
% asgard(continuity2,'timestep_method','CN')
%
% with adaptivity
% asgard(continuity2,'timestep_method','CN','adapt',true)

pde.CFL = 0.1;

soln_x = @(x,p,t)  cos(pi*x);
soln_y = @(y,p,t)  sin(2*pi*y);
soln_t = @(t)      sin(2*t);

%% Setup the dimensions
%
% Here we setup a 2D problem (x,y)

dim_x.domainMin = -1;
dim_x.domainMax = +1;
dim_x.init_cond_fn = @(x,p,t) soln_x(x,p,t) * soln_t(t);

dim_y.domainMin = -2;
dim_y.domainMax = +2;
dim_y.init_cond_fn = @(y,p,t) soln_y(y,p,t);

%%
% Add dimensions to the pde object
% Note that the order of the dimensions must be consistent with this across
% the remainder of this PDE.

pde.dimensions = {dim_x,dim_y};
num_dims = numel(pde.dimensions);

%% Setup the terms of the PDE
%
% Here we have 2 terms, having only nDims=2 (x,y) operators.

%%
% -df/dx which is 
%
% d/dx g1(x) f(x,y)          [grad,g1(x)=-1, BCL=P, BCR=P]

g1 = @(x,p,t,dat) x*0-1;
pterm1  = GRAD(num_dims,g1,0,'P','P');
term1_x = TERM_1D({pterm1});
term1   = TERM_ND(num_dims,{term1_x,[]});

%%
% -df/dy which is
%
% d/dy g1(y) f(x,y)          [grad,g1(y)=-1, BCL=P, BCR=P]

g1 = @(y,p,t,dat) y*0-1;
pterm1  = GRAD(num_dims,g1,0,'P','P');
term2_y = TERM_1D({pterm1});
term2   = TERM_ND(num_dims,{[],term2_y});

%%
% Add terms to the pde object

pde.terms = {term1,term2};

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
    @(x,p,t) cos(pi*x),   ...   % s1x
    @(y,p,t) sin(2*pi*y), ...   % s1y
    @(t)   2*cos(2*t)   ...   % s1t
    };

%%
% Source term 2
source2 = { ...
    @(x,p,t)  cos(pi*x),    ...   % s2x
    @(y,p,t)  cos(2*pi*y),  ...   % s2y
    @(t)    2*pi*sin(2*t) ...   % s2t
    };

%%
% Source term 3
source3 = { ...
    @(x,p,t)  sin(pi*x),   ...  % s3x
    @(y,p,t)  sin(2*pi*y), ...  % s3y
    @(t)    -pi*sin(2*t) ...  % s3t
    };

%%
% Add sources to the pde data structure
pde.sources = {source1,source2,source3};

%% Define the analytic solution (optional).
% This requires nDims+time function handles.

pde.analytic_solutions_1D = { ...
    soln_x,    ... % a_x
    soln_y,    ... % a_y
    soln_t     ... % a_t
    };

%%
% function to set dt

    function dt=set_dt(pde,CFL)
        
        Lmax = pde.dimensions{1}.domainMax;
        Lmin = pde.dimensions{1}.domainMin;
        LevX = pde.dimensions{1}.lev;
        dt = (Lmax-Lmin)/2^LevX*CFL;
    end

pde.set_dt = @set_dt;

end


