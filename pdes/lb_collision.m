function pde = lb_collision(opts)
% Lenard-Bernstein Equation in 1D
%
% PDE:
% 
% df/dt = div( (v-u)f + th\grad f)
%
% Domain is [0,1]x[0,1]
% Dirichlet boundary condition in advection
% Neumann boundary condition in diffusion
%
% Diffusion terms are dealt with via LDG
%
% Run with
%
% explicit
% asgard(@lb_collision,'CFL',0.01);
%
% implicit
% asgard(@lb_collision,'timestep_method','BE','dt',0.001,'num_steps',20)

%% Define the dimensions

dim_x = DIMENSION(-5,5);
dim_x.moment_dV = @(x,p,t,dat) 0*x+1;
dimensions = {dim_x};
num_dims = numel(dimensions);

%% Define the analytic solution (optional).

soln_x = @(x,p,t) p.n/sqrt(2*pi*p.th)*exp(-(x-p.u).^2/(2*p.th));
soln_t = @(t,p) 0*t+1;
soln1 = new_md_func(num_dims,{soln_x,soln_t});

solutions = {soln1};

%% Define the boundary conditions

BCL = new_md_func(num_dims,{@(x,p,t,dat) 0*x,@(t,p) 0*t+1});
BCR = new_md_func(num_dims,{@(x,p,t,dat) 0*x,@(t,p) 0*t+1});

%% Initial conditions

ic_x = @(x,p,t,d) (x > -1) - (x > 1);
ic_t = @(t,p) 0*t+1; 
ic = new_md_func(num_dims,{ic_x,ic_t});
initial_conditions = {ic};

% Define the moments of this IC
params.n = 2.011002829992263;
params.u = -4.369869711974849e-16/params.n;
params.th = 0.663189552514513/params.n+params.u.^2;

%% Construct moments

%mass moment
moment_func = new_md_func(num_dims,{@(x,p,t) 0*x+1,@(p,t)   0*t+1});
moment0 = MOMENT({moment_func});

%momentum moment
moment_func = new_md_func(num_dims,{@(x,p,t) x-p.n*p.u,@(p,t)   0*t+1});
moment1 = MOMENT({moment_func});

%energy moment
moment_func = new_md_func(num_dims,{@(x,p,t) x.^2-(p.n*p.th+p.u^2),@(p,t)   0*t+1});
moment2 = MOMENT({moment_func});

moments = {moment0,moment1,moment2};

%% Define the terms of the PDE
%
% Here we have 1 term, with each term having nDims (x and y) operators.

dV = @(x,p,t,dat) 0*x+1;

%% 
% Setup the advection div( (v-u)f )

g1 = @(x,p,t,dat) x-p.u;

pterm1 =  DIV(num_dims,g1,'',-1,'D','D',BCL,BCR,'',dV);
term1_x = SD_TERM({pterm1});

term1   = MD_TERM(num_dims,{term1_x});

%% 
% Setup the diffusion term div(theta\grad f) term

% div(theta\grad f) is the system
%
% div(sqrt(theta)q)
%
% q = sqrt(theta)\grad f

g1 = @(x,p,t,dat) x.*0+sqrt(p.th);
g2 = @(x,p,t,dat) x.*0+sqrt(p.th);

pterm1 =  DIV(num_dims,g1,'',+1,'D','D',BCL,BCR,'',dV);
pterm2 = GRAD(num_dims,g2,'',-1,'N','N','','','',dV);
term1_x = SD_TERM({pterm1,pterm2});

term2   = MD_TERM(num_dims,{term1_x});

terms = {term1, term2};

%% Define some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

%% Define sources

sources = {};

%% Define a function to set dt

    function dt=set_dt(pde,CFL)
        dims = pde.dimensions;      
        % for Diffusion equation: dt = C * dx^2
        lev = dims{1}.lev;
        dx = 1/2^lev;
        dt = CFL*dx^2;
    end

%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,[],initial_conditions,solutions,moments);

end

