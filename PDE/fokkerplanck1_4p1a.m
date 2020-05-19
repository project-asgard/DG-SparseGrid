function pde = fokkerplanck1_4p1a
% 1D test case using continuity equation, i.e., 
% df/dt == -d/dz ( (1-z^2)f )
%
% Problem is left to right convection, so we can upwind and only require
% one boundary condition, which is neumann on the left.
%
% Run with
%
% explicit
% asgard(fokkerplanck1_4p1a)

% implicit
% asgard(fokkerplanck1_4p1a,'lev',5,'deg',3,'timestep_method','CN','CFL',0.1,'num_steps',30)

pde.CFL = 0.01;

%% Setup the dimensions
% 
% Here we setup a 1D problem (x)

    function ret = phi(z,t)
        ret = tanh(atanh(z)-t);
    end
    function ret = f0(z)
        ret = z.*0+1;
    end
    function ret = soln(z,t)
        p = phi(z,t);
        t1 = 1-p.^2;
        t2 = 1-z.^2;
        t3 = f0(p);
        f = t1./t2.*t3;
        ret = f;
    end

dim_z.name = 'z';
dim_z.domainMin = -1;
dim_z.domainMax = +1;
dim_z.init_cond_fn = @(z,p,t) soln(z,0);

%%
% Add dimensions to the pde object
% Note that the order of the dimensions must be consistent with this across
% the remainder of this PDE.

pde.dimensions = {dim_z};
num_dims = numel(pde.dimensions);

%% Setup the terms of the PDE
%
% Here we have 1 term1, having only nDims=1 (x) operators.

%% 
%  -d/dz ( (1-z^2)*f )

g1 = @(z,p,t,dat) -1.*(1-z.^2);
pterm1  = GRAD(num_dims,g1,-1,'N','N');
term1_x = TERM_1D({pterm1});
term1   = TERM_ND(num_dims,{term1_x});

%%
% Add terms to the pde object

pde.terms = {term1};

%% Construct some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

pde.params = params;

%% Add an arbitrary number of sources to the RHS of the PDE
% Each source term must have nDims + 1

%%
% Add sources to the pde data structure
pde.sources = {};

%% Define the analytic solution (optional).
% This requires nDims+time function handles.

pde.analytic_solutions_1D = { ...
    @(z,p,t) soln(z,t), ...
    @(t,p) 1 
    };

%%
% Function to set time step
    function dt=set_dt(pde,CFL)
        
        dims = pde.dimensions;
        xRange = dims{1}.domainMax-dims{1}.domainMin;
        lev = dims{1}.lev;
        dx = xRange/2^lev;
        dt = CFL * dx;
        
    end
pde.set_dt = @set_dt;

end

%% Define the various input functions specified above. 


