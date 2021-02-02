function pde = advection2_spherical(opts)
% 1D test case using continuity equation, i.e., 
%
% df/dt == -div(au) = -1/r^2 d/dr( r^2 a_r u ) + 
%                           1/(r*sin(th)) d/d\th(sin(th) a_th u) = 0
% where a = [a_r,a_th] = [2,1]
%
% analytic solution 1/sqrt(r)*exp(-3*th)*cos(-.5r+t)/sin(th)
%
% Run with
%
% asgard(@advection2_spherical,'timestep_method','BE','dt',0.001,'num_steps',20,'case',1)
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
dim_x = DIMENSION(.5,2*pi);
dim_y = DIMENSION(pi/6,5*pi/6);

% dim_x = DIMENSION(0,1);

switch opts.case_
    case 1
        dim_x.jacobian = @(x,p,t) x.*0+1;
    case 2
        dim_x.jacobian = @(x,p,t) x.*0+1;
end
dimensions = {dim_x,dim_y};
num_dims = numel(dimensions);

%% Define the analytic solution (optional).

a_x = @(x,p,t) 1./sqrt(x).*cos(-.5*x+t);
a_y = @(y,p,t) exp(-3*y)./sin(y);
sol1 = new_md_func(num_dims,{a_x,a_y});
solutions = {sol1};

%% Initial conditions

initial_conditions = {sol1};

%% Define boundary conditions

BCL = sol1;
BCR = sol1;

%% Define PDE terms

switch opts.case_
    case 1     
        LHS_terms = {};
    case 2
        g1 = @(x,p,t,dat) x.^2;
        pterm1 = MASS(g1);
        g2 = @(y,p,t,dat) sin(y);
        pterm2 = MASS(g2);
        LHS_term_x = SD_TERM({pterm1});
        LHS_term_y = SD_TERM({pterm2});
        LHS_term   = MD_TERM(num_dims,{LHS_term_x,LHS_term_y});
        LHS_terms = {LHS_term};
end

% -r^2*sin(th)*a_r*u*dv/dr
g1 = @(x,p,t,dat) -x.^2*2;
g2 = @(y,p,t,dat) sin(y);
pterm1_x = GRAD(num_dims,g1,-1,'D','N',BCL);
pterm1_y = MASS(g2);
term1_x = SD_TERM({pterm1_x});
term1_y = SD_TERM({pterm1_y});
term1 = MD_TERM(num_dims,{term1_x,term1_y});

% -r*sin(th)*a_th*u*dv/dth
g1 = @(x,p,t,dat) -x;
g2 = @(y,p,t,dat) sin(y);
pterm2_x = MASS(g1);
pterm2_y = GRAD(num_dims,g2,-1,'D','N',BCL);
term2_x = SD_TERM({pterm2_x});
term2_y = SD_TERM({pterm2_y});
term2 = MD_TERM(num_dims,{term2_x,term2_y});

terms = {term1,term2};

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
