function pde = advection1_spherical(opts)
% 1D test case using continuity equation, i.e.,
%
% r^2 * df/dt == -d/dr(r^2 * v * r * f(r,t))
%
% so this is just pure advection in spherical coordinates with flux=v*r
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
% asgard(@advection1_spherical,'timestep_method','BE','dt',0.001,'num_steps',20,'case',1)
%
% case = 1 % the standard approach, i.e., no jacobian application
%               (working but minor difference near p=0)
% case = 2 % apply p^2 jacobian to LHS, volume, trace and BC evaluations using g_funcs.
%               (does not seem to work?)
% case = 3 % apply p^2 jacobian to LHS, volume, trace and BC evaluations using dim.jacobian.
%               (does not seem to work?)

%% Define paramaters

params.v = 1;

%% Define dimensions

dim_x = DIMENSION(0,2*pi);

switch opts.case_
    case 1 % no jacobian applied
        dim_x.jacobian = @(x,p,t,d) x.*0+1;
    case 2 % jacobian applied withing g_funcs
        dim_x.jacobian = @(x,p,t,d) x.*0+1;
    case 3 % jacobian applied within dim.jacobian
        dim_x.jacobian = @(x,p,t,d) x.^2;
    case 4 % jacobian applied within dim.jacobian
        dim_x.jacobian = @(x,p,t,d) x.^2;
end
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

%% Define LHS PDE terms (only on the d/dt)
switch opts.case_
    case 1 % no jacobian applied
        g1 = @(x,p,t,dat) x.^2;
    case 2 % jacobian defined within g funcs (weak form before coordinate choice)
        g1 = @(x,p,t,dat) x.^2;
    case 3 % jacobian defined within dim.jacobian (weak form before coordinate choice)
        g1 = @(x,p,t,dat) x.*0+1;
    case 4 % weak form AFTER coordinate choice
        g1 = @(x,p,t,dat) x.^2;
end
pterm1 = MASS(g1);
LHS_term_x = SD_TERM({pterm1});
LHS_term   = MD_TERM(num_dims,{LHS_term_x});
LHS_terms = {LHS_term}; % still have to include this even if it is I to allow for the jacobian application

%% Define RHS PDE terms (all other non d/dt terms)
switch opts.case_
    case 1 % no jacobian applied
        % -d/dr (r^2 * v * r * f)
        g1 = @(x,p,t,dat) -x.^2 .* p.v .* x;
    case 2 % jacobian defined within g funcs
        % -int(g*J dr) = -int(v*r * r^2 dr)
        g1 = @(x,p,t,dat) -p.v .* x .* x.^2;
    case 3 % jacobian defined within dim.jacobian
        % -int(g*J dr) = -int(v*r * J dr)
        g1 = @(x,p,t,dat) -p.v .* x;        
end

switch opts.case_
    case {1,2,3}
        pterm1 = GRAD(num_dims,g1,-1,'D','N',BCL);
        term1_x = SD_TERM({pterm1});
        term1   = MD_TERM(num_dims,{term1_x});
    case 4
        g = @(x,p,t,d) x.*0-1; 
        pterm1 = MASS(g);
        g = @(x,p,t,d) x.^2 .* x .* p.v;
        pterm2 = GRAD(num_dims,g,-1,'D','N',BCL);
        term1_x = SD_TERM({pterm1,pterm2});
        term1   = MD_TERM(num_dims,{term1_x});
end

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
pde = PDE(opts,dimensions,terms,LHS_terms,sources,params,@set_dt,[],initial_conditions,solutions);

end
