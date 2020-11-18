function pde = spacetime2(opts)
% 2D test case using continuity equation, i.e.,
%
% df/dt == 0 == - df/dx - df/dy
%
% Run with
%
% asgard(@spacetime2,'lev',4,'deg',3,'timestep_method','time_independent','num_steps',1,'adapt',true,'adapt_threshold',2e-3)

%% Define the dimensions

dim_x = DIMENSION(0,+1);
dim_y = DIMENSION(0,+1);
dimensions = {dim_x,dim_y};
num_dims = numel(dimensions);

%% Define the solutions (optional)

solutions = {};

%% Define the initial conditions

ic1 = new_md_func(num_dims); % all zero IC
initial_conditions = {ic1};

%% Define the boundary conditions

BCL = new_md_func(num_dims,{ ...
    @(x,p,t) sin(pi*x), ...
    @(y,p,t) y.*0+1, ...
    @(t,p) 1
    });

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

g1 = @(y,p,t,dat) y*0-1;
pterm1  = GRAD(num_dims,g1,0,'D','N',BCL);
term2_y = SD_TERM({pterm1});
term2   = MD_TERM(num_dims,{[],term2_y});

terms = {term1,term2};

%% Define some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

sources = {};

%% Define function to set dt

    function dt=set_dt(pde,CFL)       
        Lmax = pde.dimensions{1}.max;
        Lmin = pde.dimensions{1}.min;
        LevX = pde.dimensions{1}.lev;
        dt = (Lmax-Lmin)/2^LevX*CFL;
    end

%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,[],initial_conditions,solutions);

end


