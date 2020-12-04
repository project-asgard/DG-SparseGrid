function pde = projecti_diff1(opts)
% Example PDE using the 1D Diffusion Equation. This example PDE is
% time dependent (although not all the terms are time dependent). This
% implies the need for an initial condition.
% PDE:
%
% df/dt = nu * d^2 f/dx^2
%
% Domain is [0,2]
%
%
% d^2 f / dx^2 becomes
%
% dq/dx with free (neumann BC)
%
% and
%
% q=df/dx with Dirichlet BCs specified by analytic solution
%
% Run with
%
% explicit
% asgard(@projecti_diff1);
%
% implicit
% asgard(@projecti_diff1,'timestep_method','BE','dt',50,'num_steps',10,'lev',3,'deg',4,'case',1);


%% Define the dimensions
switch opts.case_
    case 3
        dim_x = DIMENSION(0,3);
    case 4 
        dim_x = DIMENSION(-1,2);
    case 5
        dim_x = DIMENSION(0,3);
    case 6 
        dim_x = DIMENSION(-1,2);
    otherwise
        dim_x = DIMENSION(0,2);        
end
dimensions = {dim_x};
num_dims = numel(dimensions);

%% Define the solution (optional)

k = pi/2;
nu = 0.01;
switch opts.case_
    case 1
        soln_x = @(x,p,t) sin(k*x);
    otherwise
        soln_x = @(x,p,t) cos(k*x);
end
soln_t = @(t,p) exp(-nu*k^2*t);
soln1 = new_md_func(num_dims,{soln_x,soln_t});
solutions = {soln1};

%% Define the initial condition

ic1 = soln1;
initial_conditions = {ic1};

%% Define the terms of the PDE
%
% Here we have 1 term, with each term having nDims (x and y) operators.

%%
% Setup the d^2_dx^2 term

% term1
%
% eq1 :  df/dt   == d/dx g1(x) q(x,y)   [grad,g1(x)=1, BCL=N, BCR=D]
% eq2 :   q(x,y) == d/dx g2(x) f(x,y)   [grad,g2(x)=1, BCL=D, BCR=D]
%
% coeff_mat = mat1 * mat2

g1 = @(x,p,t,dat) x.*0+nu;
g2 = @(x,p,t,dat) x.*0+1;

BCL = soln1;
BCR = soln1;

switch opts.case_
    case 1
        pterm1 = GRAD(num_dims,g1,0,'N','N'); % Q = k*cos(k*x) => Q(x) = k, Q'(x) = -k^2*sin(k*x) == 0 
        pterm2 = GRAD(num_dims,g2,0,'D','D'); % f = sin(k*x) => Q(x)=0
    case 2
        pterm1 = GRAD(num_dims,g1,+1,'D','D');
        pterm2 = GRAD(num_dims,g2,-1,'N','N');
    case 3
        pterm1 = GRAD(num_dims,g1,+1,'D','N'); % if you switch the signs of the fluxes this will fail on the right
        pterm2 = GRAD(num_dims,g2,-1,'N','D');
    case 4
        pterm1 = GRAD(num_dims,g1,-1,'N','D'); % if you switch the signs of the fluxes this will pass
        pterm2 = GRAD(num_dims,g2,+1,'D','N');
    case 5 
        pterm1 = GRAD(num_dims,g1,0,'D','N');
        pterm2 = GRAD(num_dims,g2,0,'N','D');
    case 6
        pterm1 = GRAD(num_dims,g1,0,'N','D');
        pterm2 = GRAD(num_dims,g2,0,'D','N');
end

term1_x = SD_TERM({pterm1,pterm2});
term1   = MD_TERM(num_dims,{term1_x});

terms = {term1};

%% Define some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

%% Define sources

sources = {};

    function dt=set_dt(pde,CFL)
        dims = pde.dimensions;
        % for Diffusion equation: dt = C * dx^;
        lev = dims{1}.lev;
        dx = 1/2^lev;
        dt = CFL*dx^2;
    end

%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,[],initial_conditions,solutions);

end


