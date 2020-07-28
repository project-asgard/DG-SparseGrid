function pde = fokkerplanck1_4p1b
% Problem 4.1b from the RE paper.  
% df/dt == -d/dz ( (1-z^2)f )
%
% Run with
%
% explicit
% asgard(fokkerplanck1_4p1b)
%
% implicit
% asgard(fokkerplanck1_4p1b,'timestep_method','CN','num_steps',30,'lev',4)
%
% with adaptivity
% asgard(fokkerplanck1_4p1b,'timestep_method','CN','num_steps',300,'lev',5,'deg',3,'adapt',true,'CFL',0.1)

% pde.max_lev = 8;

sig = 0.1;


%% Setup the dimensions
% 
% Here we setup a 1D problem (x)

    function ret = phi(z,t)
        ret = tanh(atanh(z)-t);
    end
    function ret = f0(z)
        ret = exp(-z.^2/sig^2);
    end
    function ret = soln(z,t)
        p = phi(z,t);
        t1 = 1-p.^2;
        t2 = 1-z.^2;
        t3 = f0(p);
        ret = t1./t2.*t3;
    end

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
% -d/dz ( (1-z^2)*f )

g1 = @(z,p,t,dat) -(1-z.^2);
pterm1  = GRAD(num_dims,g1,-1,'D','D');
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
        x_max = pde.dimensions{1}.domainMax;
        x_min = pde.dimensions{1}.domainMin;
        lev   = pde.dimensions{1}.lev;
        dt    = (x_max - x_min) / 2^lev * CFL;
    end

pde.set_dt = @set_dt;

end


