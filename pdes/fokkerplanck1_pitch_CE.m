function pde = fokkerplanck1_pitch_CE(opts)
% Problem 4.4 from the RE paper - evolution of the pitch angle dependence
% of f in the presence of electric field acceleration and collisions 
% 
% df/dt == -E d/dz((1-z^2) f) + C d/dz((1-z^2) df/dz)
%
% Run with
%
% explicit 
% asgard(@fokkerplanck1_pitch_CE,'CFL',0.01,'case',3)
%
% implicit 
% asgard(@fokkerplanck1_pitch_CE,'timestep_method','BE','num_steps',50,'dt',0.05,'lev',4,'deg',4,'case',1)

sig = 0.1;
E = 4.0;
C = 1.0;

    function ret = phi(z,t)
        ret = z.*exp(-t) ./ sqrt(1-(exp(-2*t)-1).*(z.^2));
    end
    function ret = f0(z)
        
        caseNumber = opts.case_;
        shift = 0.36;
        scaling = 0.177245;
        
        switch caseNumber
            case 1
                f = exp(-z.^2/sig^2)/scaling;
            case 2
                f = exp(-(z-shift).^2/sig^2)/scaling;
            case 3
                f = exp(-(z+shift).^2/sig^2)/scaling;
            case 4
                f = (exp(-(z-shift).^2/sig^2) + exp(-(z+shift).^2/sig^2))/(2*scaling);
        end    
        ret = f;
    end
    function ret = soln(z)
        A = E/C;
        ret = A / (2*sinh(A)) * exp(A*z);
    end

%% Define the dimensions

dim_z = DIMENSION(-1,+1);
dim_z.moment_dV = @(x,p,t,dat) 0*x+1;
dimensions = {dim_z};
num_dims = numel(dimensions);

%% Define the analytic solution (optional)

soln_z = @(z,p,t) soln(z);
soln1 = new_md_func(num_dims,{soln_z});
solutions = {soln1};

%% Define the initial conditions

ic_z = @(z,p,t) f0(z);
ic1 = new_md_func(num_dims,{ic_z});
initial_conditions = {ic1};

%% Define the terms of the PDE
%
% Here we have 1 term1, having only nDims=1 (x) operators.

%% 
% -E d/dz((1-z^2) f)
dV_z = @(z,p,t,dat) sqrt(1-z.^2);

g1 = @(z,p,t,dat) -E.*sqrt(1-z.^2);
pterm1  = DIV(num_dims,g1,'',-1,'D','D','','','',dV_z);
termE_z = SD_TERM({pterm1});
termE   = MD_TERM(num_dims,{termE_z});

%% 
% +C * d/dz( (1-z^2) df/dz )

g1 = @(z,p,t,dat) 0*z+1;
pterm1  =  DIV(num_dims,g1,'',-1,'D','D','','','',dV_z);
pterm2  = GRAD(num_dims,g1,'',+1,'N','N','','','',dV_z);
termC_z = SD_TERM({pterm1,pterm2});
termC   = MD_TERM(num_dims,{termC_z});

terms = {termE,termC};

%% Define some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

%% Define sources

sources = {};

%% Define function to set time step

    function dt=set_dt(pde,CFL)
        
        Lmax = pde.dimensions{1}.max;
        Lmin = pde.dimensions{1}.min;
        LevX = pde.dimensions{1}.lev;
        dt = (Lmax-Lmin)/2^LevX*CFL;
    end

%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,[],initial_conditions,solutions);

end

