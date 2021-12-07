function pde = continuity1(opts)
% 1D test case using continuity equation, i.e., 
%
% df/dt == -df/dx
%
% Run with
%
% explicit
% asgard(@continuity1,'lev',4,'deg',3)
%
% implicit
% asgard(@continuity1,'timestep_method','CN')

%% Define the dimensions

dim_x = DIMENSION(-1,+1);
dim_x.moment_dV = @(x,p,t,dat) 0*x+1;
dimensions = {dim_x};
num_dims = numel(dimensions);

%% Define the analytic solution (optional).
% This requires nDims+time function handles.

a_x = @(x,p,t) cos(2*pi*x);
a_t = @(t,p) sin(t);
sol1 = new_md_func(num_dims,{a_x,a_t});
solutions = {sol1};

%% Define the initial conditions

initial_conditions = {sol1};

%% Define the terms of the PDE
%
% Here we have 1 term1, having only nDims=1 (x) operators.

%% 
% -df/dx

dV = @(x,p,t,dat) 0*x+1;

g1 = @(x,p,t,dat) x.*0-1;
pterm1 = DIV(num_dims,g1,'',-1,'P','P','','','',dV);

term1_x = SD_TERM({pterm1});
term1   = MD_TERM(num_dims,{term1_x});

%%
% Add terms to the pde object

terms = {term1};

%% Define some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

%% Define sources

%%
% Source 1
s1x = @(x,p,t) cos(2*pi*x);
s1t = @(t,p) cos(t);
source1 = new_md_func(num_dims,{s1x,s1t});

%%
% Source 2
s2x = @(x,p,t) sin(2*pi*x);
s2t = @(t,p) -2*pi*sin(t);
source2 = new_md_func(num_dims,{s2x,s2t});

s3x = @(x,p,t) cos(1/2*pi*x);
s3t = @(t,p) sin(2*pi*0.1*t);
source3 = new_md_func(num_dims,{s3x,s3t});

sources = {source1,source2};

%% Define function to set time step

    function dt=set_dt(pde,CFL)
        
        dim = pde.dimensions{1};
        lev = dim.lev;
        xMax = dim.max;
        xMin = dim.min;
        xRange = xMax-xMin;
        dx = xRange/(2^lev);
        dt = CFL*dx;
    end

%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,[],initial_conditions,solutions);

end


