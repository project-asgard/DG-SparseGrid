function pde = mirror1_pitch_div(opts)
% PDE using the 1D Diffusion Equation for a case involving pitcha angle to axis \
% of the magnetic field line in a mirror configuration. This PDE is
% time dependent (although not all the terms are time dependent). The test
% particle is a hydrogen ion colliding with a background electron
% distribution that has a Maxwellian form. 
% PDE:
% 
% df/dt  == -cos(z)*f + 1/sin(z) d/dz ( sin(z) nu_D/2 df/dz ) 
%
% nu_D is the deflection frequency
% Domain is [0,pi]
% Homogeneous Neumann boundary condition 
% Code will be added to equation involving velocity dimension
%
% Diffusion term is dealt with via LDG, i.e., splitting into two first
% order equations
%
% Run with
%
% implicit
% asgard(@mirror1_pitch_div,'timestep_method','BE','dt',5e-7,'num_steps',20,'lev',3,'deg',3, 'calculate_mass', true, 'normalize_by_mass', true);
%

    function ret = psi(x)
        
        phi = erf(x);
        dphi_dx = 2./sqrt(pi) * exp(-x.^2);

        ret = 1./(2*x.^2) .* (phi - x.*dphi_dx);
    end

%Background Parameters
boltz = 1.380*10^-23; %Boltzmann constant in Joules/Kelvin
n_b = 10^19; %background density in SI units (particles/m.^3)
temp_b = 116050; %background temperature in Kelvin
% z_b = 1; %atomic number of background specie
m_b = 9.109*10^-31; %background mass in kg 
v_b = sqrt(2*boltz*temp_b/m_b); %background velocity in m/s
% eps_o = 8.85*10^-12; %permittivity of free space in Farad/m
% 
% %Target Specie Parameters
% z_a = 1;
% e = 1.602*10^-19; %charge in Coulombs
% ln_Delt = 10; %Coulomb logarithm
% m_a = 1.6726*10^-27; %target mass in kg
% coeff_ab = (e^2/(m_a*eps_o))^2; %Coefficient accounting for Coluomb force
nu_D = 10^4; %deflection frequency in s^-1

%Initial parameters for target specie
 n_o = 0.5*n_b; %initial number density for specie at specific velocity
 v_o = 0.5*v_b; %initial known velocity

%% Define the dimensions

dim_th = DIMENSION(0,pi);
dV_th = @(x,p,t,d) sin(x);
dim_th.moment_dV = dV_th;

dimensions = {dim_th};
num_dims = numel(dimensions);

%% Define the analytic solution (optional)

soln_z = @(z,p,t) n_o.*cos(z);
soln_t = @(t,p) exp(-nu_D*t);
soln1 = new_md_func(num_dims,{soln_z,soln_t});
solutions = {soln1};

%% Define the initial conditions

ic1 = soln1;
initial_conditions = {ic1};

%% Define the boundary conditions

bcl_x = @(x,p,t) x.*0;
bcl_t = soln_t;
BCL = new_md_func(num_dims,{bcl_x,bcl_t});

BCR = soln1;

%% Define the terms of the PDE
%
% Here we have 1 term1, having only nDims=1 (x) operators.

%% 

%term V2 == -cos(z)*f
%term V2 == g1(v)*f [mass, g1(v) = -cos(z), BC = N/A]
g1 = @(z,p,t,dat) -cos(z);
pterm1 = MASS(g1,'','',dV_th);
termM_z = SD_TERM({pterm1});
termM = MD_TERM(num_dims,{termM_z});

% termC is done using SLDG defining C(v) = sqrt(nu_D(v)/2)
%
% eq1 :  df/dt == div( q)   [pterm1: div (g(v)=C(v),+1, BCL=D, BCR=D]
% eq2 :      q ==  grad(f)  [pterm2: grad(g(v)=C(v),-1, BCL=N, BCR=N]


D = @(v,p) v.*0 + sqrt(nu_D/2);
g5 = @(v,p,t,dat) D(v,p);
pterm1 = DIV (num_dims,g5,'',+1,'D','D','','','',dV_th);
pterm2 = GRAD(num_dims,g5,'',-1,'N','N','','','',dV_th);
termC_z = SD_TERM({pterm1,pterm2});
termC   = MD_TERM(num_dims,{termC_z});

terms = {termM,termC};

%% Define some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

%% Define sources

sources = {};

    function dt=set_dt(pde,CFL)      
        dims = pde.dimensions;       
        % for Diffusion equation: dt = C * dx^2       
        lev = dims{1}.lev;
        dx = 1/2^lev;
        dt = CFL*dx^2;        
    end

%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,[],initial_conditions,solutions);

end