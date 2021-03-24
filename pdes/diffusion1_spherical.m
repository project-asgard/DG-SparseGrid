function pde = diffusion1_spherical(opts)
% Example of 1D diffusion in a non-cartesian coordinate system. This will
% test the non-unit mass martix inversion within the LDG method. 
% 
% p^2 * df/dt == -d/dp(p^2 * df/dp)
%
% Run with
%
% asgard(@diffusion1_spherical,'timestep_method','BE');
%
% Notes
%

%% Define the dimensions
dim_x = DIMENSION(0,1);
dim_x.jacobian = @(x,p,t,d) x.^2;
dimensions = {dim_x};
num_dims = numel(dimensions);

%% Define the analytic solution
soln_x = @(x,p,t) exp(-x.^2);
soln1 = new_md_func(num_dims,{soln_x});
solutions = {soln1};

%% Initial conditions
initial_conditions = {soln1};

%% Define the LHS terms (the d/dt)

gL = @(x,p,t,dat) x.^2;
ptermL = MASS(gL);
termL_p = SD_TERM({ptermL});
termL   = MD_TERM(num_dims,{termL_p});
LHSterms = {termL};

%% Define the RHS terms (everything except d/dt)

% LHS = -d/dp(p^2 * df/dp)
%
% eq1 :  LHS == d/dp(g1*w)   [grad,g1=-p^2, BCL=N, BCR=D]
% eq2 :    w == d/dp(g2*f)   [grad,g2=1, BCL=D, BCR=D]
%
% coeff_mat = mat1 * mat2

g1 = @(x,p,t,dat) -x.^2;
g2 = @(x,p,t,dat) x.*0+1;

pterm1 = GRAD(num_dims,g1,+1,'D','N');
pterm2 = GRAD(num_dims,g2,-1,'N','D');

term1_x = SD_TERM({pterm1,pterm2});
term1   = MD_TERM(num_dims,{term1_x});

terms = {term1};

%% Define some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

%% Define sources

s1x = @(x,p,t) exp(-x.^2).*(-6+4.*x.^2).*x.^2.*x.^2; % note the extra x^2 here is the jacobian
s1t = @(t,p) t.*0+1;
source1 = {s1x, s1t};

sources = {source1};

    function dt=set_dt(pde,CFL)
        
        dims = pde.dimensions;
        
        % for Diffusion equation: dt = C * dx^2
        
        lev = dims{1}.lev;
        dx = 1/2^lev;
        dt = CFL*dx^2;
        
    end

%% Construct PDE

pde = PDE(opts,dimensions,terms,LHSterms,sources,params,@set_dt,[],initial_conditions,solutions);

end
