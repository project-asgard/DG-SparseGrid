function pde = continuity1(opts)
% 1D test case using continuity equation, i.e., 
%
% df/dt == -df/dx
%
% Run with
%
% explicit
% asgard(@continuity1)
% asgard(@continuity1,'lev',4,'deg',3)
%
% implicit
% asgard(@continuity1,'timestep_method','CN')
% asgard(@continuity1,'timestep_method','CN','CFL',0.1)

%% Define the dimensions

dim_x = DIMENSION(-1,+1);
dimensions = {dim_x};
num_dims = numel(dimensions);

%% Define the initial conditions

ic_x = @(x,p,t) x.*0;
ic1 = new_md_func(num_dims,{ic_x});
initial_conditions = {ic1};

%% Define the terms of the PDE
%
% Here we have 1 term1, having only nDims=1 (x) operators.

%% 
% -df/dx

g1 = @(x,p,t,dat) x.*0-1;
pterm1 = GRAD(num_dims,g1,0,'P','P');

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
source1 = {s1x,s1t};

%%
% Source 2
s2x = @(x,p,t) sin(2*pi*x);
s2t = @(t,p) -2*pi*sin(t);
source2 = {s2x,s2t};

sources = {source1,source2};

%% Define the analytic solution (optional).
% This requires nDims+time function handles.

a_x = @(x,p,t) cos(2*pi*x);
a_t = @(t,p) sin(t);
sol1 = new_md_func(num_dims,{a_x,a_t});
solutions = {sol1};

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


