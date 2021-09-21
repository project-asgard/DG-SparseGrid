function pde = mirror1_space_div(opts)
% 1D mirror case along magnetic axis direction, i.e., 
%
% df/dt f B(s) == -cos(pi/4)*df/ds B(s) - cos(pi/4)*dB/ds/B(s) f B(s)
%
% Run with
%
% explicit
% asgard(@mirror1_space_div, 'calculate_mass', false, 'normalize_by_mass', false)

% implicit
% asgard(@mirror1_space_div,'timestep_method','BE', 'calculate_mass', true, 'normalize_by_mass', true)

params = mirror_parameters();

params.B_func = @(s) exp(s); %magnetic field as a function of spatial coordinate
params.dB_ds = @(s) params.B_func(s); %derivative of magnetic field

params.vel_test = 1; %test value for velocity in coefficient
params.pitch_test = pi/4; %test value for pitch angle in coefficient
decay_coeff = @(v,z) v.*cos(z);

%% Define the dimensions

dim_s = DIMENSION(-3,3);
dV_s = @(x,p,t,d) params.B_func(x); 
dim_s.moment_dV = dV_s;
dimensions = {dim_s};
num_dims = numel(dimensions);

%% Define the analytic solution

sol_s = @(s,p,t) exp(s);
sol_t = @(t,p) exp(-2*decay_coeff(params.vel_test,params.pitch_test).*t);
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
% -B(s)*d(vtest*cos(ztest)f)/ds
g1 = @(s,p,t,dat) s.*0 - decay_coeff(p.vel_test,p.pitch_test)*(p.pitch_test > pi/2);
pterm1 = DIV(num_dims,g1,'',-1,'D','D','','','', dV_s);
%pterm1 = GRAD(num_dims,g1,-1,'N','N');
term1a_s = SD_TERM({pterm1});
term1a   = MD_TERM(num_dims,{term1a_s});

g1 = @(s,p,t,dat) s.*0 - decay_coeff(p.vel_test,p.pitch_test)*(p.pitch_test < pi/2);
pterm1 = DIV(num_dims,g1,'',+1,'D','D','','','', dV_s);
%pterm1 = GRAD(num_dims,g1,-1,'N','N');
term1b_s = SD_TERM({pterm1});
term1b   = MD_TERM(num_dims,{term1b_s});

%%
% Add terms to the pde object

%% Mass term
% B(s)*(vcosz/B)dB/ds f
g2 = @(s,p,t,dat) s.*0 - decay_coeff(p.vel_test,p.pitch_test).*(p.dB_ds(s)./p.B_func(s)).*(p.pitch_test < pi/2).*(p.B_func(s) > 0);
pterm1   = MASS(g2,'','',dV_s);
termB1a = SD_TERM({pterm1});

term2a = MD_TERM(num_dims,{termB1a});

g2 = @(s,p,t,dat) s.*0 - decay_coeff(p.vel_test,p.pitch_test).*(p.dB_ds(s)./p.B_func(s)).*(p.pitch_test > pi/2).*(p.B_func(s) < 0);
pterm1   = MASS(g2,'','',dV_s);
termB1b = SD_TERM({pterm1});

term2b = MD_TERM(num_dims,{termB1b});

terms = {term1a,term1b,term2a,term2b};

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
