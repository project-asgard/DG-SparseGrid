function pde = fokkerplanck1_momentum_E_div(opts)
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
% p is a spherical coordinate (p,th,ph), so the div and grad look like 
%
% div[] = 1/p^2 * d/dp * p^2[], grad[] = d/dp[]
%
% and the volument_element dV = p^2
%
% setting E=2 and z=0.5, and expanding the d/dz term via product rule, then setting d/dz=0,
% we get for the momentum direction ...
%
% termE1 == -1/p^2 d/dp p^2 f(p)           [pterm1: div(g(p)=-1,+1, BCL=?, BCR=?]
% termE2 == 2/p f(p)                       [pterm1: mass(g(p)=2/p]
%
% Run with
%
% asgard(@fokkerplanck1_momentum_E_div,'timestep_method','matrix_exponential','dt',0.01,'num_steps',10,'case',2)
%
% case = 1 % flat initial condition
% case = 2 % maxwellian initial condition

params = fokkerplanck_parameters(opts);

%% Setup the dimensions 

dV_p = @(x,p,t,d) x.^2;
dim_p = DIMENSION(0.1,+10);
dim_p.volume_element = dV_p;
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
        ic_p = f0_p;
        ic1 = new_md_func(num_dims,{ic_p});
        initial_conditions = {ic1};
    case 2
        initial_conditions = {soln1};
end

%% Define the terms of the PDE

%% -div(flux_E) == termE1

% termE1 == -1/p^2 (d/dp p^2 f(p)) [div, g1(p) = -1,   BCL=D,BCR=D]

g1     = @(x,p,t,d) x.*0-1; 
pterm1 = DIV (num_dims,g1,'',-1,'D','N',soln1,soln1,'',dV_p);

termE1_p = SD_TERM({pterm1});
termE1   = MD_TERM(num_dims,{termE1_p});

% termE2 == +2/p f
%        == g1(p) f(p)       [mass, g1(p) = +2/p,  BC N/A]

g1       = @(x,p,t,d) +2./x;
pterm1   = MASS(g1,'','',dV_p);
termE2_p = SD_TERM({pterm1});
termE2   = MD_TERM(num_dims,{termE2_p});

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


