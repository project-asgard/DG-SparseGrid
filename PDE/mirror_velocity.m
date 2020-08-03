function pde = mirror_velocity
% One-dimensional magnetic mirror from the FP paper - evolution of the velocity dependence
% of f in the presence of electric field acceleration and collisions 
% 
% df/dt == -(Z_a E/m_a)d/dv(f) + 1/v^2 (d/dv(flux_v))
%
% flux_v == v^3[(m_a/(m_a + m_b))nu_s f) + 0.5*nu_par*v*d/dv(f)]
%
% Run with
%
% asgard(mirror_velocity,'timestep_method','BE', 'dt', 1e-8);

pde.CFL = 0.01;
%Background Parameters
k_b = 1.380*10^-23; %Boltzmann constant in Joules/Kelvin
n_b = 10^19; %background density in SI units (particles/m.^3)
T_b = 116050; %background temperature in Kelvin
z_b = 1; %atomic number of background specie
m_b = 9.109*10^-31; %background mass in kg 
v_b = sqrt(2*k_b*T_b/m_b); %background velocity in m/s
eps_o = 8.85*10^-12; %permittivity of free space in Farad/m

%Target Specie Parameters
z_a = 1;
e = 1.602*10^-19; %charge in Coulombs
ln_Delt = 10; %Coulomb logarithm
m_a = 1.6726*10^-27; %target mass in kg
L_ab = (e^2/(m_a*eps_o))^2; %Coefficient accounting for Coluomb force
nu_s = 10^6; %Slowing down frequency in s^-1
nu_par = 10^5; %parallel diffusion frequency


%E = 1.0; %parallel Electric field

%% Setup the dimensions
% 
function y = dd1(n)
% Our default value is 0
y = 0; 

% The function is 1 only if the input is 0
if n == 0
    y = 1;
end

end
dim_v.domainMin = 0.1;
dim_v.domainMax = 2*v_b;
dim_v.init_cond_fn = @(v,p,t) dd1(v-v_b);

%%
% Add dimensions to the pde object
% Note that the order of the dimensions must be consistent with this across
% the remainder of this PDE.

pde.dimensions = {dim_v};
num_dims = numel(pde.dimensions);

%% Setup the terms of the PDE

%% 
% -E*Z_a/m_a d/dz(f)

%g1 = @(v,p,t,dat) -E.*Z_a/m_a;
%pterm1  = GRAD(num_dims,g1,-1,'N','N');
%termE_v = TERM_1D({pterm1});
%termE   = TERM_ND(num_dims,{termE_v});

%% 
% d/dv(v^3(m_a/(m_a + m_b))nu_s f))

g1 = @(v,p,t,dat) v.^3*m_a*nu_s/(m_a + m_b);
pterm1  = GRAD(num_dims,g1,-1,'N','N');
termV_s = TERM_1D({pterm1});
termV1   = TERM_ND(num_dims,{termV_s});

%%
% d/dv(v^4*0.5*nu_par*d/dv(f))
% d/dv(g(v)r(v))
% r(v) = d/dv(f)

g1 = @(v,p,t,dat) v.^4.*0.5*nu_par;
g2 = @(v,p,t,dat) v.*0 + 1;
pterm1 = GRAD(num_dims,g1,-1,'D','D');
pterm2 = GRAD(num_dims,g2,+1,'N','N');
termV_par = TERM_1D({pterm1,pterm2});
termV2   = TERM_ND(num_dims,{termV_par});

%%
% Add terms to the pde object

pde.terms = {termV1,termV2};

%% Construct some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

pde.params = params;

%% Add an arbitrary number of sources to the RHS of the PDE
% Each source term must have nDims + 1

%%
% Add sources to the pde data structure
pde.sources = {};

%% Define the analytic solution (optional).
% This requires nDims+time function handles.

pde.analytic_solutions_1D = { ...
    @(z,p,t) z.*0 + 1, ...
    @(t,p) t.*0 + 1 
    };

%%
% Function to set time step
    function dt=set_dt(pde,CFL)
        
        Lmax = pde.dimensions{1}.domainMax;
        LevX = pde.dimensions{1}.lev;
        dt = Lmax/2^LevX*CFL;
    end

pde.set_dt = @set_dt;

end