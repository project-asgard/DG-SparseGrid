function pde = fokkerplanck2_E(opts)
% Combining momentum and pitch angle dynamics for the E term
%
% d/dt f(p,z) == -div(flux_E)
%
% flux_E is the flux due to E accleration
%
% -div(flux_E) == termE1 + termE2
%
% termE1 == -E*z*f(z) * 1/p^2 (d/dp p^2 f(p))
% termE2 == -E/p*f(p) * d/dz (1-z^2) f(z)
%
% setting E=2 and z=0.5, and expanding the d/dz term, then setting d/dz=0,
% we get for the momentum direction ...
%
% termE1 = -1/p^2 d/dp p^2 f(p) + 2/p f(p)
%
% Run with
%
% asgard(@fokkerplanck1_momentum_E,'timestep_method','matrix_exponential','dt',0.01,'num_steps',10,'case',2)
%
% case = 1 % flat initial condition
% case = 2 % maxwellian initial condition

params = fokkerplanck_parameters(opts);

%% Setup the dimensions 

dim_p = DIMENSION(0.0,+1);
dim_p.jacobian = @(x,p,t) x.^2;
dimensions = {dim_p};
num_dims = numel(dimensions);

%% Define the analytic solution (optional)

switch opts.case_
    case 1
        soln_p = @(x,p,t) x.*0 + 1;
    case 2
        soln_p = @(x,p,t) exp(-(x-t).^2);
end

soln1 = new_md_func(num_dims,{soln_p});
solutions = {soln1};

%% Define the initial conditions

switch opts.case_
    case 1
        f0_p = @(x,p,t) x.*0+1;
    case 2
        f0_p = @(x,p,t) exp(-x.^2);
end

ic_p = f0_p;
ic1 = new_md_func(num_dims,{ic_p});
initial_conditions = {ic1};

%% Define the terms of the PDE

%% -div(flux_E) == termE1 + termE2

% termE1 == -1/p^2 (d/dp p^2 f(p))
%        == q(p)
%   q(p) == g2(p) u(p)       [mass, g2(p) = 1/p^2, BC N/A]
%   u(p) == d/dp g3(p) f(p)  [grad, g3(p) = p^2,   BCL=N,BCR=D]

g2 = @(x,p,t,dat) -min(1./x.^2,1e12);
g3 = @(x,p,t,dat) x.^2;
pterm1   = MASS(g2); 
pterm2   = GRAD(num_dims,g3,+1,'N','N');% Lin's Setting (DLG-why is this central flux?)
termE1_p = SD_TERM({pterm1,pterm2});
termE1 = MD_TERM(num_dims,{termE1_p});

% termE2 == +2/p f
%        == g1(p) f(p)       [mass, g1(p) = +2/p,  BC N/A]

g1 = @(x,p,t,dat) +2./x;
pterm1   = MASS(g1);
termE2_p = SD_TERM({pterm1});
termE2= MD_TERM(num_dims,{termE2_p});

terms = {termE1,termE2};


%% Define sources

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


