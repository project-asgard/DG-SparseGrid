function pde = diffusion1_div(opts)
% Example PDE using the 1D Diffusion Equation. This example PDE is
% time dependent (although not all the terms are time dependent). This
% implies the need for an initial condition. 
% PDE:
% 
% df/dt = d^2 f/dx^2
%
% Domain is [0,1]
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
% asgard(@diffusion1,'CFL',0.01);
%
% implicit
% asgard(@diffusion1,'timestep_method','BE');
%
% Notes
%
% Louis: With small nu, stability is maintained for longer times
%   asgard(@diffusion1,'lev',3,'deg',4,'timestep_method','BE', 'dt',0.05,'num_steps',20)
% DLG : it seems like the k and nu are mixed up / together here? i.e., the
% equation assumes nu = 1 since it is not present in g1 or g2, but nu is
% specified as != 1 below?

%% Define the dimensions
% 
% Here we setup a 1D problem (x,y)
nu = pi/2; %coefficient set to be very small to allow for stability % DLG - WTF is this?
soln_x = @(x,p,t) cos(nu*x);
%soln_x = @(x,p,t) sin(nu*x);
soln_t = @(t,p) exp(-2*nu^2*t);

BCFunc = @(x,p,t) soln_x(x,p,t);
BCFunc_t = @(t,p) soln_t(t,p);

% Domain is (a,b)

% The function is defined for the plane
% x = a and x = b
BCL_fList = { ...
    @(x,p,t) BCFunc(x,p,t), ... % replace x by a
    @(t,p) BCFunc_t(t,p)
    };

BCR_fList = { ...
    @(x,p,t) BCFunc(x,p,t), ... % replace x by b
    @(t,p) BCFunc_t(t,p)
    };

dim_x = DIMENSION(0,1);
dim_x.moment_dV = @(x,p,t,dat) 0*x+1;
dimensions = {dim_x};
num_dims = numel(dimensions);

%% Initial conditions

ic1 = new_md_func(num_dims,{soln_x,soln_t});
initial_conditions = {ic1};

%% Define the terms of the PDE
%
% Here we have 1 term, with each term having nDims (x and y) operators.

%% 
% Setup the d^2_dx^2 term

% term1
%
% eq1 :  df/dt   == d/dx g1(x) q(x,y)   [div ,g1(x)=1, BCL=N, BCR=D]
% eq2 :   q(x,y) == d/dx g2(x) f(x,y)   [grad,g2(x)=1, BCL=D, BCR=D]
%
% coeff_mat = mat1 * mat2

dV = @(x,p,t,dat) 0*x+1;

g1 = @(x,p,t,dat) x.*0+1;
g2 = @(x,p,t,dat) x.*0+1;

pterm1 =   DIV(num_dims,g1,'',-1,'N','N','','','',dV);
pterm2 =  GRAD(num_dims,g2,'',+1,'D','D',BCL_fList,BCR_fList,'',dV);

term1_x = SD_TERM({pterm1,pterm2});
term1   = MD_TERM(num_dims,{term1_x});

%% Hack for interior penalty
pen = 0;
g3 = @(x,p,t,dat) 0.*x+pen;
g4 = @(x,p,t,dat) 0.*x-pen;

pterm1 = DIV(num_dims,g3,'',-1,'D','D','','','',dV);
term2_x = SD_TERM({pterm1});
term2 = MD_TERM(num_dims,{term2_x});

pterm1 = DIV(num_dims,g4,'', 0,'D','D','','','',dV);
term3_x = SD_TERM({pterm1});
term3 = MD_TERM(num_dims,{term3_x});

%%% With penalty
terms = {term1,term2,term3};

%%% Without penalty
%terms = {term1};

%% Define some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

%% Define sources

s1x = @(x,p,t) -nu^2*cos(nu*x);
s1t = @(t,p) exp(-2*nu^2*t);
source1 = new_md_func(num_dims,{s1x, s1t});

%sources = {};
sources = {source1};

%% Define the analytic solution (optional).
% This requires nDims+time function handles.

soln1 = new_md_func(num_dims,{soln_x,soln_t});
solutions = {soln1};

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
