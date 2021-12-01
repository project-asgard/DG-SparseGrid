function pde = advection1_div(opts)
% 1D test case using continuity equation, i.e., 
%
% df/dt == -2*df/dx - 2*sin(x)
%
% Run with
%
% explicit
% asgard(@advection1_div)
% asgard(@advection1,'lev',4,'deg',3)
%
% implicit
% asgard(@advection1_div,'timestep_method','CN')
% asgard(@advection1_div,'timestep_method','CN','CFL',0.01)


%% Define dimensions

dim_x = DIMENSION(0,pi);
dim_x.moment_dV = @(x,p,t,d) 0*x+1;
dimensions = {dim_x};
num_dims = numel(dimensions);

%% Initial conditions
ic_x = @(x,p,t) cos(x);
ic1 = new_md_func(num_dims,{ic_x});
initial_conditions = {ic1};

%% Define boundary conditions

% Dirichlet on the right, f(0) = 1

BCFunc_Left  = @(x) x.*0 +1;
BCFunc_Right = @(x) x.*0 -1;

BCL_fList = { ...
    @(x,p,t) BCFunc_Left(x), ... % replace x by b
    @(t,p) t.*0 + 1
    };

BCR_fList = { ...
    @(x,p,t) BCFunc_Right(x), ... % replace x by b
    @(t,p) t.*0 + 1
    };


%% Define PDE terms
% Here we have 1 term1, having only nDims=1 (x) operators.

dV = @(x,p,t,dat) 0*x+1;
 
% -2*df/dx
g1 = @(x,p,t,dat) x.*0-2;
pterm1 = GRAD(num_dims,g1,'',-1,'D','N', BCL_fList, '','',dV);

term1_x = SD_TERM({pterm1});
term1   = MD_TERM(num_dims,{term1_x});

terms = {term1};


%% Define paramaters
%  These might be used within the various functions.

params.parameter1 = 0;
params.parameter2 = 1;


%% Define sources

% Source 1
s1x = @(x,p,t) -2.*sin(x);
s1t = @(t,p) t.*0 + 1;
source1 = {s1x,s1t};

sources = {source1};


%% Define the analytic solution (optional).

a_x = @(x,p,t) cos(x);
sol1 = new_md_func(num_dims,{a_x});
solutions = {sol1};

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
