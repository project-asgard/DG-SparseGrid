function pde = fokkerplanck1_4p5
% Problem 4.5 from the RE paper - evolution of the pitch angle dependence
% of f in the presence of electric field acceleration and collisions and
% radiation damping
% 
% df/dt == -E d/dz((1-z^2) f) + C d/dz((1-z^2) df/dz) - R d/dz(z(1-z^2) f)
%
% Run with
%
% asgard(fokkerplanck1_4p5,'implicit',true,'lev',3,'CFL',0.5,'num_steps',50)

pde.CFL = 0.01;
sig = 0.1;
E = 2.0;
C = 1.0;
R = 2.0;

%% Setup the dimensions
% 
% Here we setup a 1D problem (x)

    function ret = f0(z)
        
        caseNumber = 1;
        
        switch caseNumber
            case 1
                f = exp(-z.^2/sig^2);
        end    
        ret = f;
    end
    function ret = soln(z,t)
        A = E/C;
        B = R/C;
        Q = .03;
        ret = Q * exp(A*z + (B/2)*z.^2);
    end

dim_z.domainMin = -1;
dim_z.domainMax = +1;
dim_z.init_cond_fn = @(z,p,t) f0(z);

%%
% Add dimensions to the pde object
% Note that the order of the dimensions must be consistent with this across
% the remainder of this PDE.

pde.dimensions = {dim_z};
num_dims = numel(pde.dimensions);

%% Setup the terms of the PDE

%% 
% -E d/dz((1-z^2) f)

g1 = @(z,p,t,dat) -E.*(1-z.^2);
pterm1  = GRAD(num_dims,g1,-1,'D','D');
termE_z = TERM_1D({pterm1});
termE   = TERM_ND(num_dims,{termE_z});

%% 
% +C * d/dz( (1-z^2) df/dz )

g1 = @(z,p,t,dat) 1-z.^2;
g2 = @(z,p,t,dat) z.*0+1;
pterm1  = GRAD(num_dims,g1,-1,'D','D');
pterm2  = GRAD(num_dims,g2,+1,'N','N');
termC_z = TERM_1D({pterm1,pterm2});
termC   = TERM_ND(num_dims,{termC_z});

%%
% - R d/dz(z(1-z^2) f)

g1 = @(z,p,t,dat) -R * z.*(1-z.^2);
pterm1  = GRAD(num_dims,g1,-1,'D','D');
termR_z = TERM_1D({pterm1});
termR   = TERM_ND(num_dims,{termR_z});

%%
% Add terms to the pde object

pde.terms = {termE,termC,termR};

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
        
        Lmax = pde.dimensions{1}.domainMax;
        LevX = pde.dimensions{1}.lev;
        dt = Lmax/2^LevX*CFL;
    end

pde.set_dt = @set_dt;

end

