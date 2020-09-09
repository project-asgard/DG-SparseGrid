function pde = mirror_space1
% 1D mirror case along magnetic axis direction, i.e., 
%
% df/dt == -vcosz*df/ds -(vcosz/B)dB/ds f
%
% Run with
%
% explicit
% asgard(mirror_space1)

% implicit
% asgard(advection1,'mirror_space1','BE')


%% Setup the dimensions
% 
% Here we setup a 1D problem (s)

B_func = @(s) exp(s); %magnetic field as a function of spatial coordinate
dB_ds = @(s) B_func(s); %derivative of magnetic field

v_test = 500; %test value for velocity in coefficient
z_test = pi/2 - 1e-6; %test value for pitch angle in coefficient
A = v_test*cos(z_test);

%% Functions that define the analytic solution
a_s = @(s) exp(s);
a_t = @(t) exp(-2*A*t);

dim_s.domainMin = 0;
dim_s.domainMax = 5;
dim_s.init_cond_fn = @(s,p,t) a_s(s);

%%

% Add dimensions to the pde object
% Note that the order of the dimensions must be consistent with this across
% the remainder of this PDE.

pde.dimensions = {dim_s};
num_dims = numel(pde.dimensions);

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


%% Setup the terms of the PDE
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

%%Mass term
%(vcosz/B)dB/ds f
g2 = @(s,p,t,dat) s.*0 - A.*dB_ds(s)./B_func(s);
pterm1   = MASS(g2);
termB1 = TERM_1D({pterm1});

term2 = TERM_ND(num_dims,{termB1});

pde.terms = {term1,term2};

%% Construct some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

pde.params = params;

%% Add an arbitrary number of sources to the RHS of the PDE
% Each source term must have nDims + 1 (which here is 2+1 / (x,v) + time) functions describing the
% variation of each source term with each dimension and time.
% Here we define 3 source terms.

%%
% Source 1
s1x = @(x,p,t) -2.*sin(x);
s1t = @(t,p) t.*0 + 1;
source1 = {s1x,s1t};

%%
% Source 2
%s2x = @(x,p,t) sin(2*pi*x);
%s2t = @(t,p) -2*pi*sin(t);
%source2 = {s2x,s2t};

%%
% Add sources to the pde data structure
pde.sources = {};

%% Define the analytic solution (optional).
% This requires nDims+time function handles.

pde.analytic_solutions_1D = { ...
    @(s,p,t) a_s(s), ...
    @(t,p) a_t(t)
    };

%%
% Function to set time step

    function dt=set_dt(pde,CFL)
        
        dim = pde.dimensions{1};
        lev = dim.lev;
        xMax = dim.domainMax;
        xMin = dim.domainMin;
        xRange = xMax-xMin;
        dx = xRange/(2^lev);
        dt = CFL*dx;
    end

pde.set_dt = @set_dt;

end