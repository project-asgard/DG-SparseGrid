function pde = diffusion2(opts)
% Example PDE using the 2D (1x-1y) Heat Equation. This example PDE is
% time dependent (although not all the terms are time dependent). This
% implies the need for an initial condition. 
% PDE:
% 
% df/dt = d^2 f/dx^2 + d^2 f/dy^2
%
% Domain is [0,1]x[0,1]
% Dirichlet boundary condition 
%
% Diffusion terms are dealt with via LDG, i.e., splitting into two first
% order equations:
%
% d^2 f / dx^2 becomes
%
% dq/dx with free (homogeneous Neumann BC)
%
% and
%
% q=df/dx with Dirichlet BCs specified by analytic solution
%
% Run with
%
% explicit
% asgard(@diffusion2,'CFL',0.01);
%
% implicit
% asgard(@diffusion2,'timestep_method','BE','dt',0.001,'num_steps',20)

%% Define the dimensions

dim_y = DIMENSION(0,1);
dim_x = DIMENSION(0,1);
dimensions = {dim_x,dim_y};
num_dims = numel(dimensions);

%% Define the analytic solution (optional).

soln_x = @(x,p,t) cos(pi*x);
soln_y = @(y,p,t) cos(pi*y);
soln_t = @(t,p) exp(-2*pi^2*t);
soln1 = new_md_func(num_dims,{soln_x,soln_y,soln_t});

solutions = {soln1};

%% Define the boundary conditions

BCL = soln1;
BCR = soln1;

%% Initial conditions

initial_conditions = {soln1};

%% Define the terms of the PDE
%
% Here we have 1 term, with each term having nDims (x and y) operators.

%% 
% Setup the d^2_dx^2 term

% term1
%
% eq1 :  df/dt   == d/dx g1(x) q(x,y)   [grad,g1(x)=2, BCL=N, BCR=N]
% eq2 :   q(x,y) == d/dx g2(x) f(x,y)   [grad,g2(x)=1, BCL=D, BCR=D]
%
% coeff_mat = mat1 * mat2

g1 = @(x,p,t,dat) x.*0+1;
g2 = @(x,p,t,dat) x.*0+1;

pterm1 = GRAD(num_dims,g1,+1,'N','N');
pterm2 = GRAD(num_dims,g2,-1,'D','D',BCL,BCR);

term1_x = SD_TERM({pterm1,pterm2});
term1   = MD_TERM(num_dims,{term1_x,[]});

%% 
% Setup the d^2_dy^2 term

% term2
%
% eq1 :  df/dt   == d/dy g1(y) q(x,y)   [grad,g1(y)=2, BCL=N, BCR=N]
% eq2 :   q(x,y) == d/dy g2(y) f(x,y)   [grad,g2(y)=1, BCL=D, BCR=D]
%
% coeff_mat = mat1 * mat2

g1 = @(y,p,t,dat) y.*0+1;
g2 = @(y,p,t,dat) y.*0+1;

pterm1 = GRAD(num_dims,g1,+1,'N','N');
pterm2 = GRAD(num_dims,g2,-1,'D','D',BCL,BCR);

term2_y = SD_TERM({pterm1,pterm2});
term2   = MD_TERM(num_dims,{[],term2_y});

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

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,[],initial_conditions,solutions);

end

