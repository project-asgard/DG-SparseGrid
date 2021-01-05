function pde = mirror_space1(opts)
% 1D mirror case along magnetic axis direction, i.e., 
%
% df/dt == -vcosz*df/ds -(vcosz/B)dB/ds f
%
% Run with
%
% explicit
% asgard(@mirror_space1)

% implicit
% asgard(@mirror_space1,'timestep_method','BE')

B_func = @(s) exp(s); %magnetic field as a function of spatial coordinate
dB_ds = @(s) B_func(s); %derivative of magnetic field

vel_test = 500; %test value for velocity in coefficient
pitch_test = pi/2 - 1e-6; %test value for pitch angle in coefficient
decay_coeff = vel_test*cos(pitch_test);

%% Define the dimensions

dim_s = DIMENSION(0,5);
dimensions = {dim_s};
num_dims = numel(dimensions);

%% Define the analytic solution

sol_s = @(s,p,t) exp(s);
sol_t = @(t,p) exp(-2*decay_coeff*t);
soln1 = new_md_func(num_dims,{sol_s,sol_t});
solutions = {soln1};

%% Define the initial conditions

ic1 = new_md_func(num_dims,{sol_s});
initial_conditions = {ic1};

%% Define the boundary conditions

BCL = soln1;
BCR = soln1;

%% Define the terms of the PDE
%
% We have an advection and mass term

%% Advection  term
% -vcosz*df/ds
g1 = @(s,p,t,dat) s.*0 - decay_coeff;
pterm1 = GRAD(num_dims,g1,-1,'D','D', BCL, BCR);
%pterm1 = GRAD(num_dims,g1,-1,'N','N');
term1_s = SD_TERM({pterm1});
term1   = MD_TERM(num_dims,{term1_s});

%%
% Add terms to the pde object

%% Mass term
% (vcosz/B)dB/ds f
g2 = @(s,p,t,dat) s.*0 - decay_coeff.*dB_ds(s)./B_func(s);
pterm1   = MASS(g2);
termB1 = SD_TERM({pterm1});

term2 = MD_TERM(num_dims,{termB1});

terms = {term1,term2};

%% Define some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

%% Define sources

sources = {};

%% Define function to set time step

    function dt=set_dt(pde,CFL)      
        dim = pde.dimensions{1};
        lev = dim.lev;
        xMax = dim.max;
        xMin = dim.min;
        xRange = xMax-xMin;
        dx = xRange/(2^lev);
        dt = CFL*dx;
    end

%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,[],initial_conditions,solutions);

end
