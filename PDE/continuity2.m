function pde = continuity2(opts)
% 2D test case using continuity equation, i.e.,
%
% df/dt == -df/dx - df/dy
%
% Run with
%
% explicit
% asgard(@continuity2,'lev',3,'deg',3,'CFL',0.1)
%
% implicit
% asgard(@continuity2,'timestep_method','CN')
%
% with adaptivity
% asgard(@continuity2,'timestep_method','CN','adapt',true)

% pde.CFL = 0.1;

soln_x = @(x,p,t)  cos(pi*x);
soln_y = @(y,p,t)  sin(2*pi*y);
soln_t = @(t)      sin(2*t);

%% Define the dimensions
%
% Here we setup a 2D problem (x,y)

dim_x = DIMENSION(-1,+1);
dim_x.init_cond_fn = @(x,p,t) soln_x(x,p,t) * soln_t(t);
dim_x.jacobian = @(x,p,t) x.*0 + 1;

dim_y = DIMENSION(-2,+2);
dim_y.init_cond_fn = @(y,p,t) soln_y(y,p,t);
dim_y.jacobian = @(y,p,t) y.*0 + 1;

dimensions = {dim_x,dim_y};
num_dims = numel(dimensions);

%% Define the terms of the PDE
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

terms = {term1,term2};

%% Define some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;


%% Define sources

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

sources = {source1,source2,source3};

%% Define the analytic solution (optional).
% This requires nDims+time function handles.

analytic_solutions_1D = { ...
    soln_x,    ... % a_x
    soln_y,    ... % a_y
    soln_t     ... % a_t
    };

%% Define function to set dt

    function dt=set_dt(pde,CFL)
        
        Lmax = pde.dimensions{1}.max;
        Lmin = pde.dimensions{1}.min;
        LevX = pde.dimensions{1}.lev;
        dt = (Lmax-Lmin)/2^LevX*CFL;
    end

%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,analytic_solutions_1D);

end


