function pde = advection1_bidirectional_flow(opts)
% 1D test case for where the sign of advection changes within the domain 
%
% df/dt == -d/dx (x*f)
%
% Run with
%
% implicit
% asgard(@advection1_bidirectional_flow,'timestep_method','BE','dt',0.01,'num_steps',50,'deg',3)


%% Define dimensions

dim_x = DIMENSION(-pi,+pi);
dimensions = {dim_x};
num_dims = numel(dimensions);

%% Define the analytic solution (optional).

a_x = @(x,p,t) exp(-t) .* cos(exp(-t).*x);
sol1 = new_md_func(num_dims,{a_x});
solutions = {sol1};

%% Initial conditions

ic1 = sol1;
initial_conditions = {ic1};

%% Define boundary conditions

BCL = sol1;
BCR = sol1;

%% Define PDE terms
 
% -d/dx (x*f)
g1 = @(x,p,t,dat) -x;
pterm1 = GRAD(num_dims,g1,-1,'D','D', BCL, BCR);

term1_x = SD_TERM({pterm1});
term1   = MD_TERM(num_dims,{term1_x});

terms = {term1};

%% Define paramaters

params.parameter1 = 0;
params.parameter2 = 1;

%% Define sources

sources = {};

%% Define function to calculate time step

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
