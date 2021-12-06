function pde = advection1_reverse_flow_div(opts)
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
dim_x.moment_dV = @(x,p,t,dat) 0*x+1;

dimensions = {dim_x};
num_dims = numel(dimensions);

%% Define initial conditions

ic_x = @(x,p,t) cos(x);
ic1 = new_md_func(num_dims,{ic_x,[]});

initial_conditions = {ic1};

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

dV = @(x,p,t,dat) 0*x+1;

g1 = @(x,p,t,dat) x.*0 + 2;
pterm =  DIV(num_dims,g1,'',-1,'N','D', '', BCR_fList,'',dV);
term_x = SD_TERM({pterm});
term1 = MD_TERM(num_dims, {term_x});

terms = {term1};

%% Define parameters
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

%% Define sources

s1x = @(x,p,t) 2*sin(x);
source1 = new_md_func(num_dims,{s1x,[]});

sources = {source1};


%% Define the analytic solution (optional).
% This requires nDims+time function handles.
a_x = @(x,p,t) cos(x);
soln1 = new_md_func(num_dims,{a_x,[]});

solutions = {soln1};

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

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,[],initial_conditions,solutions);

end
