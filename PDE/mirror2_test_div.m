function pde = mirror2_test_div(opts)

% Two-dimensional test of spatial and pitch-angle variation in magnetic
% mirror where f coefficients change signs throughout the domain
% 
% df/dt ==  1/(sin(th)) d/dth (sin(th) sin(th)/2B dB/ds f_a)
%           - 1/B d (B cos(th)f)/ds
% Run with
%
% asgard(@mirror2_space_pitch_div,'timestep_method','BE','case',3,'dt',1e-8,'lev',3,'deg', 3,'num_steps',20)

dim_th = DIMENSION(0,pi);
dV_th = @(x,p,t,d) sin(x);
dim_th.moment_dV = dV_th;

dimensions = {dim_th};
num_dims = numel(dimensions);

params = mirror_parameters();

%% Define the initial conditions

params.init_cond_z = @(z,p,t,dat) cos(z);
ic1 = new_md_func(num_dims,{params.init_cond_z});
initial_conditions = {ic1};

%% Define the boundary conditions

BCL = new_md_func(num_dims,{...
    @(x,p,t,dat) cos(x), ...
    @(t,p) 0*t+1});

BCR = new_md_func(num_dims,{...
    @(x,p,t,dat) cos(x), ...
    @(t,p) 0*t+1});

%% Define the terms of the PDE

%% Advection  term

%termZa is a done combining mass ad div defining C(th) = sin(th) and using
%the function mag(s) = dB/ds/(2*B). This applies to the domain where 
% mag(s) > 0 and therefore we use upwinding
%eq1: df/dt == div(mag(s) f)       [pterm1: div(g(th)=C(th),-1, BCL=N,BCR=N)]

dV_th = @(x,p,t,d) sin(x);
%dV_th = @(x,p,t,d) 0*x+1;

C = @(z) sin(z);
%C = @(z) 1./sin(z);
%C = @(z) 0*z+1;
g2 = @(z,p,t,dat) C(z);
pterm2 = DIV(num_dims,g2,'',-1,'N','D','','','',dV_th);
termZa_z = SD_TERM({pterm2});

termZa = MD_TERM(num_dims,{termZa_z});

g2 = @(z,p,t,dat) -cos(z);
pterm2 = MASS(g2,'','',dV_th);
termZa_mass = SD_TERM({pterm2});
term_mass = MD_TERM(num_dims,{termZa_mass});

terms = {termZa,term_mass};
%terms = {termZa};



%% Define the analytic solution (optional)


soln1_z = @(z,p,t,d) cos(z);
solution = new_md_func(num_dims,{soln1_z});
solutions = {solution};


%% Define sources
%source1_v = @(x,p,t,d) p.f0_v(x);
%source1_z = @(x,p,t,d) (sin(x)).^2 - (cos(x)).^2;
%source1 = new_md_func(num_dims,{source1_z,[]});
%sources = {source1};
sources = {};

%% Define function to set time step
    function dt=set_dt(pde,CFL)      
        Lmax = pde.dimensions{1}.max;
        LevX = pde.dimensions{1}.lev;
        dt = Lmax/2^LevX*CFL;
    end



%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,[],initial_conditions,solutions);
end