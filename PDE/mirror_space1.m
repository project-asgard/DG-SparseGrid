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

v_test = 500; %test value for velocity in coefficient
z_test = pi/2 - 1e-6; %test value for pitch angle in coefficient
A = v_test*cos(z_test);

%% Define the analytic solution

a_s = @(s) exp(s);
a_t = @(t) exp(-2*A*t);

%% Define the dimensions

dim_s = DIMENSION(0,5);
dim_s.init_cond_fn = @(s,p,t) a_s(s);

dimensions = {dim_s};
num_dims = numel(dimensions);

% Setup boundary conditions of the solution
% Dirichlet on the right, f(0) = 1

%BCFunc_t = @(t) soln_t(t);

BCL_fList = { ...
    @(s,p,t) a_s(s), ... % replace x by b
    @(t,p) a_t(t)
    };

BCR_fList = { ...
    @(s,p,t) a_s(s), ... % replace x by b
    @(t,p) a_t(t)
    };


%% Define the terms of the PDE
%
% We have an advection and mass term

%% Advection  term
% -vcosz*df/ds
g1 = @(s,p,t,dat) s.*0 - A;
pterm1 = GRAD(num_dims,g1,-1,'D','D', BCL_fList, BCR_fList);

term1_s = TERM_1D({pterm1});
term1   = TERM_ND(num_dims,{term1_s});

%%
% Add terms to the pde object

%% Mass term
% (vcosz/B)dB/ds f
g2 = @(s,p,t,dat) s.*0 - A.*dB_ds(s)./B_func(s);
pterm1   = MASS(g2);
termB1 = TERM_1D({pterm1});

term2 = TERM_ND(num_dims,{termB1});

terms = {term1,term2};

%% Define some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

%% Define sources

sources = {};

%% Define the analytic solution (optional).
% This requires nDims+time function handles.

analytic_solutions_1D = { ...
    @(s,p,t) a_s(s), ...
    @(t,p) a_t(t)
    };

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

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,analytic_solutions_1D);

end
