function pde = advec_diffusion1
% Example PDE using the 1D Advection-Diffusion Equation. This example PDE is
% time dependent (although not all the terms are time dependent). This
% implies the need for an initial condition. 
% PDE:
% 
% df/dt = -df/dx + d^2 f/dx^2 + (cos(x) - sin(x))exp(-t)
%
% Domain is [0, pi/4]
% Inhomogeneous Dirichlet boundary condition 
% Code is expected to be added to the 1d Advection Equation 
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
% asgard(advec_diffusion1);
%
% implicit
% asgard(advec_diffusion1,'implicit',true);

pde.CFL = 0.01;

%% Setup the dimensions
% 
% Here we setup a 1D problem (x,y)

soln_x = @(x) cos(x) + sin(x);
soln_t = @(t) exp(-t);

BCFunc = @(x) soln_x(x);
BCFunc_t = @(t) soln_t(t);

% Domain is (a,b)

% The function is defined for the plane
% x = a and x = b
BCL_fList = { ...
    @(x,p,t) BCFunc(x), ... % replace x by a
    @(t,p) BCFunc_t(t)
    };

BCR_fList = { ...
    @(x,p,t) BCFunc(x), ... % replace x by b
    @(t,p) BCFunc_t(t)
    };
%Incorrect BC to check whether each term requires BC
%Each term DOES require correct BC
%BCIncorrect_fList= { ...
%    @(x,p,t) x.*0 + 2, ... % replace x by b
%    @(t,p) BCFunc_t(t)
%    };

dim_x.name = 'x';
dim_x.domainMin = 0;
dim_x.domainMax = pi/4;
dim_x.init_cond_fn = @(x,p,t) soln_x(x)*soln_t(t);


%%
% Add dimensions to the pde object
% Note that the order of the dimensions must be consistent with this across
% the remainder of this PDE.

pde.dimensions = {dim_x};
num_dims = numel(pde.dimensions);

%% Setup the terms of the PDE
%
% Here we have 2 terms, with each term having nDims (x) operators.

%%
% Setup the -d_dx term

g1 = @(x,p,t,dat) x.*0-1;
pterm1 = GRAD(num_dims,g1,-1,'D','D', BCL_fList, BCR_fList);

term1_x = TERM_1D({pterm1});
term1   = TERM_ND(num_dims,{term1_x});

%% 
% Setup the d^2_dx^2 term

% term1
%
% eq1 :  df/dt   == d/dx g1(x) q(x,y)   [grad,g1(x)=1, BCL=N, BCR=N]
% eq2 :   q(x,y) == d/dx g2(x) f(x,y)   [grad,g2(x)=1, BCL=D, BCR=D]
%
% coeff_mat = mat1 * mat2

g1 = @(x,p,t,dat) x.*0+1;
g2 = @(x,p,t,dat) x.*0+1;

pterm1 = GRAD(num_dims,g1,+1,'N','N');
pterm2 = GRAD(num_dims,g2,-1,'D','D', BCL_fList, BCR_fList);

term2_x = TERM_1D({pterm1,pterm2});
term2   = TERM_ND(num_dims,{term2_x});

%%
% Add terms to the pde object

 pde.terms = {term1, term2};

%% Construct some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

pde.params = params;

%%
% Sources

s1x = @(x,p,t) cos(x) - sin(x);
s1t = @(t,p) exp(-t);
source1 = {s1x, s1t};

pde.sources = {source1};

%% Define the analytic solution (optional).
% This requires nDims+time function handles.

pde.analytic_solutions_1D = { ...
    @(x,p,t) soln_x(x), ...
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
