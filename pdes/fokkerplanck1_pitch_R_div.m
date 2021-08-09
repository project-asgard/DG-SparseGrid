function pde = fokkerplanck1_pitch_R_div(opts)
% Problem 4.3 from the RE paper - radiation damping term  
% df/dt == -d/dz ( z(1-z^2)f )
%
% df/dt = -div( z*sqrt(1-z^2)f\hat{z} ) 
%
% Run with
%
% explicit
% asgard(@fokkerplanck1_pitch_R_div,'CFL',0.01)
%
% implicit
% asgard(@fokkerplanck1_pitch_R_div,'timestep_method','matrix_exponential','lev',4,'deg',4,'num_steps',20,'dt',0.01)

sig = 0.1;

%% Define the dimensions
% 
% Here we setup a 1D problem (x)

    function ret = phi(z,t)
        ret = z.*exp(-t) ./ sqrt(1+(exp(-2*t)-1).*(z.^2));
    end
    function ret = f0(z)
        
        shift = 0.36;
        
        switch opts.case_
            case 1
                f = exp(-z.^2/sig^2);
            case 2
                f = exp(-(z-shift).^2/sig^2);
            case 3
                f = exp(-(z+shift).^2/sig^2);
            case 4
                f = exp(-(z-shift).^2/sig^2) + exp(-(z+shift).^2/sig^2);
        end    
        ret = f;
    end
    function ret = soln(z,t)
        p = phi(z,t);
        t1 = p.*(1-p.^2);
        t2 = z.*(1-z.^2);
        t3 = f0(p);
        ret = t1./t2.*t3;
    end

dim_z = DIMENSION(-1,+1);
dim_z.moment_dV = @(x,p,t,dat) 0*x+1;
dimensions = {dim_z};
num_dims = numel(dimensions);

%% Define the analytic solution (optional)

soln_z = @(z,p,t) soln(z,t);
soln1 = new_md_func(num_dims,{soln_z});
solutions = {soln1};

%% Define the initial conditions

ic_z = @(z,p,t) soln(z,0);
ic1 = new_md_func(num_dims,{ic_z});
initial_conditions = {ic1};

%% Define the terms of the PDE
%
% Here we have 1 term1, having only nDims=1 (x) operators.

%% 
% -d/dz ( z(1-z^2)f )

g1 = @(z,p,t,dat) -z.*sqrt(1-z.^2);
dV = @(z,p,t,dat) sqrt(1-z.^2);
pterm1  = DIV(num_dims,g1,'',-1,'D','D','','','',dV);
term1_z = SD_TERM({pterm1});
term1   = MD_TERM(num_dims,{term1_z});

terms = {term1};

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

