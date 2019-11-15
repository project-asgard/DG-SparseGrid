function pde = advection1reverse
%1D Test case of the 1D advection equation with
% the direction of motion reversed from advection1.m.
% Code designed to test inhomogeneous Dirichlet BC in conjunction with
% advection1.m
% Equation to be solved is 
% df/dt = df/dx + source(x)
% Run with
%
% explicit
% asgard(advection1reverse)
% asgard(advection1reverse,'lev',4,'deg',3,'implicit',false)
%
% implicit
% asgard(advection1reverse,'implicit',true)
% asgard(advection1reverse,'implicit',true,'CFL',0.01)

%% Setup the dimensions
% 
% Here we setup a 1D problem (x) on [a,b]

dim_x.domainMin = -pi/2; %a
dim_x.domainMax = 0; %b
dim_x.init_cond_fn = @(x,p,t) cos(x);

%%
% Add dimensions to the pde object
% Note that the order of the dimensions must be consistent with this across
% the remainder of this PDE.

pde.dimensions = {dim_x};
num_dims = numel(pde.dimensions);

%Inhomogeneous Dirichlet condition on one side of the domain

BC_Func = @(x) x.*0 + 1;
%BCFunc_t = @(t) soln_t(t);

BCL_fList = { ... 
     @(x,p,t) x.*0, ... %replace x with a
     @(t,p)  t.*0 + 1 %boundary condition for time, usually 1
     };

BCR_fList = { ... 
     @(x,p,t) BC_Func(x), ... %replace x with b
     @(t,p)  t.*0 + 1 %boundary condition for time, usually 1
     };

%Setup the terms of the PDE
% Here the term is df/dx, which is opposite the direction in advection1.m

g1rev = @(x,p,t,dat) x.*0 + 1;
pterm = GRAD(num_dims,g1rev,-1,'D','D', BCL_fList, BCR_fList); %flux is in opposite direction of advection1.m

term_x = TERM_1D({pterm});
term1 = TERM_ND(num_dims, {term_x});

%Add terms to the pde object
pde.terms = {term1};

%% Construct some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

pde.params = params;

%Add the source term, making sure to note for the sign difference in the example equation
% Each source term must have nDims + 1 (which here is 2+1 / (x,v) + time) functions describing the
% variation of each source term with each dimension and time.

%Source
s1x = @(x,p,t) sin(x);
s1t = @(t,p) t.*0 + 1;
source1 = {s1x, s1t};

%Add sources to the pde data structure
pde.sources = {source1};


%% Define the analytic solution (optional).
% This requires nDims+time function handles.
a_x = @(x,p,t) cos(x);
a_t = @(t,p) t.*0 + 1;

pde.analytic_solutions_1D = {a_x, a_t};

%function to set time step
    function dt = set_dt(pde, CFL)
        dim = pde.dimensions{1};
        lev = dim.lev;
        xMax = dim.domainMax;
        xMin = dim.domainMin;
        xRange = xMax - xMin;
        dx = xRange/(2^lev);
        dt = CFL*dx;
    end

pde.set_dt = @set_dt;

end
