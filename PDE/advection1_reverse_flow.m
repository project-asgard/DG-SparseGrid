function pde = advection1_reverse_flow(opts)
% 1D advection problem with inhomogeneous Dirichlet BC 
% (flow direction reversed from advection1.m)
%  
% df/dt = 2*df/dx + 2*sin(x)
%
% Run with
%
% explicit
% asgard(@advection1_reverse_flow)
% asgard(@advection1_reverse_flow,'lev',4,'deg',3)
%
% implicit
% asgard(@advection1_reverse_flow,'timestep_method', 'CN')
% asgard(@advection1_reverse_flow,'timestep_method', 'CN','CFL',0.01)

%% Define dimensions

dim_x = DIMENSION(0,pi);
dim_x.init_cond_fn = @(x,p,t) cos(x);

dimensions = {dim_x};

num_dims = numel(dimensions);

%Inhomogeneous Dirichlet condition on one side of the domain

BC_Func_Left  = @(x) x.*0 +1;
BC_Func_Right = @(x) x.*0 -1;

BCL_fList = { ... 
     @(x,p,t) BC_Func_Left(x), ... %replace x with a
     @(t,p)  t.*0 + 1 %boundary condition for time, usually 1
     };

BCR_fList = { ... 
     @(x,p,t) BC_Func_Right(x), ... %replace x with b
     @(t,p)  t.*0 + 1 %boundary condition for time, usually 1
     };

%% Define PDE terms
% Here the term is df/dx, which is opposite the direction in advection1.m

g1 = @(x,p,t,dat) x.*0 + 2;
pterm = GRAD(num_dims,g1,-1,'D','D', BCL_fList, BCR_fList);

term_x = TERM_1D({pterm});
term1 = TERM_ND(num_dims, {term_x});

terms = {term1};

%% Define parameters
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

%% Define sources

s1x = @(x,p,t) 2*sin(x);
s1t = @(t,p) t.*0 + 1;
source1 = {s1x, s1t};

sources = {source1};


%% Define the analytic solution (optional).
% This requires nDims+time function handles.
a_x = @(x,p,t) cos(x);
a_t = @(t,p) t.*0 + 1;

analytic_solutions_1D = {a_x, a_t};

%% Define function to set time step
    function dt = set_dt(pde, CFL)
        dim = pde.dimensions{1};
        lev = dim.lev;
        xMax = dim.max;
        xMin = dim.min;
        xRange = xMax - xMin;
        dx = xRange/(2^lev);
        dt = CFL*dx;
    end

%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,analytic_solutions_1D);

end
