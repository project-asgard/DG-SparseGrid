function pde = advection1_spherical(opts)
% 1D test case using continuity equation, i.e., 
%
% df/dt == -1/r^2 * d/dr r^2 * v * r * f(r,t)
%
% f(r,0) == cos(r)
%
% analytic solution exp(-3*t*v) * cos(exp(-t*v)*r)
%
% Run with
%
% asgard(@advection1_spherical,'timestep_method','BE','dt',0.001,'num_steps',20,'case',1)
%
% case = 1 % the standard approach, i.e., no jacobian application 
%               (working but minor difference near p=0)
% case = 2 % apply p^2 jacobian to LHS, volume, trace and BC evaluations.
%               (does not seem to work?)

%% Define paramaters

params.v = 1;

%% Define dimensions

%dim_x = DIMENSION(0.5,2*pi);
%dim_x = DIMENSION(1000,1001);
dim_x = DIMENSION(0,2*pi);

% dim_x = DIMENSION(0,1);

switch opts.case_
    case 1
        dim_x.jacobian = @(x,p,t) x.*0+1;
    case 2
        dim_x.jacobian = @(x,p,t) x.^2;
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

%% Define PDE terms

switch opts.case_
    case 1     
        g1 = @(x,p,t,dat) x.^2;
        pterm1 = MASS(g1);
        LHS_term_x = SD_TERM({pterm1});
        LHS_term   = MD_TERM(num_dims,{LHS_term_x});
        LHS_terms = {LHS_term};
        
        pterm = GRAD(num_dims,@(x,p,t,dat) -p.v .* x.^3,-1,'D','N',BCL);
        term1_x = SD_TERM({pterm});
        term1 = MD_TERM(num_dims,{term1_x});
        
    case 2
        g1 = @(x,p,t,dat) 0*x+1;
        pterm1 = MASS(g1);
        LHS_term_x = SD_TERM({pterm1});
        LHS_term   = MD_TERM(num_dims,{LHS_term_x});
        LHS_terms = {LHS_term};
        
        pterm = GRAD(num_dims,@(x,p,t,dat) -p.v .* x,-1,'D','N',BCL);
        term1_x = SD_TERM({pterm});
        term1 = MD_TERM(num_dims,{term1_x});
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
