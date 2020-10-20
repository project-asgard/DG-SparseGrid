function pde = fokkerplanck1_4p5(opts)
% Problem 4.5 from the RE paper - evolution of the pitch angle dependence
% of f in the presence of electric field acceleration and collisions and
% radiation damping
% 
% df/dt == -E d/dz((1-z^2) f) + C d/dz((1-z^2) df/dz) - R d/dz(z(1-z^2) f)
%
% Run with
%
% asgard(@fokkerplanck1_4p5,'timestep_method','CN','lev',3,'CFL',0.5,'num_steps',50,'case',1)

sig = 0.1;
E = 2.0;
C = 1.0;
R = 2.0;

%% Define the dimensions
% 
% Here we setup a 1D problem (x)

    function ret = f0(z)
        
        caseNumber = opts.case_;
        
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

dim_z = DIMENSION(-1,+1);
dim_z.init_cond_fn = @(z,p,t) f0(z);

dimensions = {dim_z};
num_dims = numel(dimensions);

%% Define the terms of the PDE

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

terms = {termE,termC,termR};

%% Define some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;


%% Define the sources

sources = {};

%% Define the analytic solution (optional).
% This requires nDims+time function handles.

analytic_solutions_1D = { ...
    @(z,p,t) soln(z,t), ...
    @(t,p) 1 
    };

%% Define the function to set time step
    function dt=set_dt(pde,CFL)      
        Lmax = pde.dimensions{1}.max;
        LevX = pde.dimensions{1}.lev;
        dt = Lmax/2^LevX*CFL;
    end

%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,analytic_solutions_1D);

end

