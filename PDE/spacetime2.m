function pde = spacetime2(opts)
% 2D test case using continuity equation, i.e.,
%
% df/dt == 0 == - df/dx - df/dy
%
% Run with
%
% asgard(@spacetime2,'lev',4,'deg',3,'timestep_method','time_independent','num_steps',1,'adapt',true,'adapt_threshold',2e-3)

soln_x = @(x,p,t)  x.*0;
soln_y = @(y,p,t)  y.*0;
soln_t = @(t)      1;

%% Setup the dimensions
%
% Here we setup a 2D problem (x,y)

dim_x = DIMENSION(0,+1);
dim_x.init_cond_fn = @(x,p,t) x.*0;

dim_y = DIMENSION(0,+1);
dim_y.init_cond_fn = @(y,p,t) y.*0;

dimensions = {dim_x,dim_y};
num_dims = numel(dimensions);

%% Define the terms of the PDE
%
% Here we have 2 terms, having only nDims=2 (x,y) operators.

%%
% -df/dx which is 
%
% d/dx g1(x) f(x,y)          [grad,g1(x)=-1, BCL=D, BCR=N]

g1 = @(x,p,t,dat) x*0-1;
pterm1  = GRAD(num_dims,g1,0,'D','N');
term1_x = SD_TERM({pterm1});
term1   = MD_TERM(num_dims,{term1_x,[]});

%%
% -df/dy which is
%
% d/dy g1(y) f(x,y)          [grad,g1(y)=-1, BCL=D, BCR=N]
% 
BCL_fList = { ...
    @(x,p,t) sin(pi*x), ... % replace x by a
    @(y,p,t) y.*0+1, ...
    @(t,p) 1
    };

g1 = @(y,p,t,dat) y*0-1;
pterm1  = GRAD(num_dims,g1,0,'D','N',BCL_fList);
term2_y = SD_TERM({pterm1});
term2   = MD_TERM(num_dims,{[],term2_y});

terms = {term1,term2};

%% Define some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

sources = {};

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


