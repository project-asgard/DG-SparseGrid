function pde = advection1_spherical_div(opts)
% 1D test case using continuity equation, i.e.,
%
% df/dt == -1/r^2*d/dr(r^2 * v * r * f(r,t))
%
% so this is just pure advection in spherical coordinates with flux=v*r*f
%
% f(r,0) == cos(r)
%
% direction of advection is always to the right (+ve), so left BC is D,
% while right is free (N)
%
% analytic solution exp(-3*t*v) * cos(exp(-t*v)*r)
%
% Run with
%
% asgard(@advection1_spherical,'timestep_method','BE','dt',0.001,'num_steps',20)
%

%% Define paramaters

params.v = 1;

%% Define dimensions

dim_x = DIMENSION(0,2*pi);
dim_x.moment_dV = @(x,p,t,dat) x.^2;
dimensions = {dim_x};
num_dims = numel(dimensions);

%% Define the analytic solution (optional).
a_x = @(x,p,t) exp(-3.*t.*p.v) .* cos(exp(-t.*p.v).*x);
sol1 = new_md_func(num_dims,{a_x});
solutions = {sol1};

%% Initial conditions
initial_conditions = {sol1};

%% Define boundary conditions
BCL = sol1;
BCR = sol1;

%% Define PDE terms (all other non d/dt terms)

dV = @(x,p,t,dat) x.^2;

g1 = @(x,p,t,dat) -p.v*x;
pterm1 = DIV(num_dims,g1,'',-1,'D','N','','','',dV); 
        %Note left BC does not matter since jacobian is zero there.
        %Right BC is free since flow is going to the right
term1_x = SD_TERM({pterm1});
term1   = MD_TERM(num_dims,{term1_x});

terms = {term1};

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
