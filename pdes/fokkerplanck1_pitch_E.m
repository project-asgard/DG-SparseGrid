function pde = fokkerplanck1_pitch_E(opts)
% 1D test case using continuity equation, i.e., 
% df/dt == -d/dz ( (1-z^2)f )
%
% df/dt = -div( sqrt(1-z^2)f\hat{z} )
%
% Problem is left to right convection, so we can upwind and only require
% one boundary condition, which is Dirichlet on the left.
%
% When the domain is [-1,1] then the BC don't really matter (as they will
% be annihilated) 
%
% Run with
%
% explicit
% asgard(@fokkerplanck1_pitch_E,'case',1)

% implicit
% asgard(@fokkerplanck1_pitch_E,'lev',5,'deg',3,'timestep_method','CN','CFL',0.1,'num_steps',30,'case',1)
%
% case = 1; % flat initial condition
% case = 2; % gaussian initial condition

    function ret = phi(z,t)
        ret = tanh(atanh(z)-t);
    end
    function ret = f0(z)
        case_number = opts.case_;
        switch case_number
            case 1
                ret = z.*0+1;
            case 2
                sig = 0.1;
                ret = exp(-z.^2/sig^2);
        end
    end
    function ret = soln(z,t)
        p = phi(z,t);
        t1 = 1-p.^2;
        t2 = 1-z.^2;
        t3 = f0(p);
        ret = t1./t2.*t3;
    end

%% Define the dimensions

dim_z = DIMENSION(-1,+1);
dim_z.moment_dV = @(x,p,t,dat) 0*x+1;
dimensions = {dim_z};
num_dims = numel(dimensions);

%% Define the analytic solution (optional).

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
%  -d/dz ( (1-z^2)*f )

g1 = @(z,p,t,dat) -sqrt(1-z.^2);
dV_z = @(z,p,t,dat) sqrt(1-z.^2);

pterm1  = DIV(num_dims,g1,'',-1,'N','N','','','',dV_z);
term1_x = SD_TERM({pterm1});
term1   = MD_TERM(num_dims,{term1_x});

terms = {term1};

%% Define some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

%% Define source

sources = {};

%% Define function to set time step
    function dt=set_dt(pde,CFL)      
        dims = pde.dimensions;
        xRange = dims{1}.max-dims{1}.min;
        lev = dims{1}.lev;
        dx = xRange/2^lev;
        dt = CFL * dx;     
    end

%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,[],initial_conditions,solutions);

end

