function pde = fokkerplanck1_pitch_CER(opts)
% Problem 4.5 from the RE paper - evolution of the pitch angle dependence
% of f in the presence of electric field acceleration and collisions and
% radiation damping
%
% df/dt == -E d/dz((1-z^2) f) + C d/dz((1-z^2) df/dz) - R d/dz(z(1-z^2) f)
%
% df/dt == -div( (F_E + F_C + F_R)\hat{z} )
%
% where
%
% F_E = E*sqrt(1-z^2)f
%
% F_C = C*grad(f)
%
% F_R = R*z*sqrt(1-z^2)f
%
% Run with
%
% asgard(@fokkerplanck1_pitch_CER,'timestep_method','matrix_exponential','lev',3,'deg',4,'dt',1,'num_steps',50,'case',1,'normalize_by_mass',true)

sig = 0.1;
E = 2.0;
C = 1.0;
R = 2.0;

    function ret = f0(z)
        
        caseNumber = opts.case_;
        
        switch caseNumber
            case 1
                f = exp(-z.^2/sig^2);
        end
        ret = f;
    end
    function ret = soln(z,p,t)
        A = E/C;
        B = R/C;
        Q = .03;
        ret = Q * exp(A*z + (B/2)*z.^2);
        if isfield(p,'norm_fac')
            ret = p.norm_fac .* ret;
        end
    end

%% Define the dimensions

dim_z = DIMENSION(-1,+1);
dim_z.moment_dV = @(x,p,t,dat) 0*x+1;
dimensions = {dim_z};
num_dims = numel(dimensions);

%% Define the initial conditions

ic_z = @(z,p,t) f0(z);
ic1 = new_md_func(num_dims,{ic_z});
initial_conditions = {ic1};

%% Define the terms of the PDE
dV_z = @(z,p,t,dat) sqrt(1-z.^2);

%%
% -E d/dz((1-z^2) f)

g1 = @(z,p,t,dat) -E.*sqrt(1-z.^2);
pterm1  =  DIV(num_dims,g1,'',-1,'D','D','','','',dV_z);
termE_z = SD_TERM({pterm1});
termE   = MD_TERM(num_dims,{termE_z});

%%
% +C * d/dz( (1-z^2) df/dz )

g1 = @(z,p,t,dat) 0*z+1;
pterm1  =  DIV(num_dims,g1,'',-1,'D','D','','','',dV_z);
pterm2  = GRAD(num_dims,g1,'',+1,'N','N','','','',dV_z);
termC_z = SD_TERM({pterm1,pterm2});
termC   = MD_TERM(num_dims,{termC_z});

%%
% - R d/dz(z(1-z^2) f)

g1 = @(z,p,t,dat) -R * z.*sqrt(1-z.^2);
pterm1  =  DIV(num_dims,g1,'',-1,'D','D','','','',dV_z);
termR_z = SD_TERM({pterm1});
termR   = MD_TERM(num_dims,{termR_z});

terms = {termE,termC,termR};

%% Define some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

%% Define the sources

sources = {};

%% Define the analytic solution (optional).

soln_z = @(z,p,t) soln(z,p,t);
soln1 = new_md_func(num_dims,{soln_z});
solutions = {soln1};

%% Define the function to set time step
    function dt=set_dt(pde,CFL)
        Lmax = pde.dimensions{1}.max;
        LevX = pde.dimensions{1}.lev;
        dt = Lmax/2^LevX*CFL;
    end

%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,[],initial_conditions,solutions);

end

