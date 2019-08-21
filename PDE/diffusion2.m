function pde = diffusion2
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
% asgard(diffusion2);
%
% implicit
% asgard(diffusion2,'implicit',true);

pde.CFL = 0.01;

%% Setup the dimensions
% 
% Here we setup a 2D problem (x,y)

soln_x = @(x) cos(pi*x);
soln_y = @(y) cos(pi*y);
soln_t = @(t) exp(-2*pi^2*t);

BCFunc = @(x) soln_x(x);
BCFunc_t = @(t) soln_t(t);

% Domain is (a,b)x(c,d)

% The function is defined for the plane
% x = a and x = b
BCL_fList = { ...
    @(x,p,t) BCFunc(x), ... % replace x by a
    @(y,p,t) BCFunc(y), ...
    @(t,p) BCFunc_t(t)
    };

BCR_fList = { ...
    @(x,p,t) BCFunc(x), ... % replace x by b
    @(y,p,t) BCFunc(y), ...
    @(t,p) BCFunc_t(t)
    };

dim_x.name = 'x';
dim_x.domainMin = 0;
dim_x.domainMax = 1;
dim_x.init_cond_fn = @(x,p,t) soln_x(x)*soln_t(t);

% The function is defined for the plane
% y = c and y = d
BCL_fList = { ...
    @(x,p,t) BCFunc(x), ...
    @(y,p,t) BCFunc(y), ... % replace y by c
    @(t,p) BCFunc_t(t)
    };

BCR_fList = { ...
    @(x,p,t) BCFunc(x), ...
    @(y,p,t) BCFunc(y), ...  % replace y by d
    @(t,p) BCFunc_t(t)
    };

dim_y.name = 'y';
dim_y.domainMin = 0;
dim_y.domainMax = 1;
dim_y.init_cond_fn = @(y,p,t) soln_y(y)*soln_t(t);

%%
% Add dimensions to the pde object
% Note that the order of the dimensions must be consistent with this across
% the remainder of this PDE.

pde.dimensions = {dim_x, dim_y};
num_dims = numel(pde.dimensions);

%% Setup the terms of the PDE
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
pterm2 = GRAD(num_dims,g2,-1,'D','D',BCL_fList,BCR_fList);

term1_x = TERM_1D({pterm1,pterm2});
term1   = TERM_ND(num_dims,{term1_x,[]});

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
pterm2 = GRAD(num_dims,g2,-1,'D','D',BCL_fList,BCR_fList);

term2_y = TERM_1D({pterm1,pterm2});
term2   = TERM_ND(num_dims,{[],term2_y});

%%
% Add terms to the pde object

 pde.terms = {term1, term2};

%% Construct some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

pde.params = params;

%%
% Add sources to the pde data structure
pde.sources = {};

%% Define the analytic solution (optional).
% This requires nDims+time function handles.

pde.analytic_solutions_1D = { ...
    @(x,p,t) soln_x(x), ...
    @(y,p,t) soln_y(y), ... 
    @(t,p) soln_t(t) 
    };

    function dt=set_dt(pde,CFL)
        
        dims = pde.dimensions;
        
        % for Diffusion equation: dt = C * dx^2
        
        lev = dims{1}.lev;
        dx = 1/2^lev;
        dt = CFL*dx^2;
        
    end

pde.set_dt = @set_dt;

end

%%
% Function to set time step

