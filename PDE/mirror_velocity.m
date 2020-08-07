function pde = mirror_velocity
% One-dimensional magnetic mirror from the FP paper - evolution of the velocity dependence
% of f in the presence of electric field acceleration and collisions 
% 
% df/dt == -(Z_a E/m_a)d/dv(f) + 1/v^2 (d/dv(flux_v))
%a
% flux_v == v^3[(m_a/(m_a + m_b))nu_s f) + 0.5*nu_par*v*d/dv(f)]
%
% Run with
%
% asgard(mirror_velocity,'timestep_method','BE', 'dt', 1e-9, 'num_steps', 100, 'lev', 6, 'deg', 7)

pde.CFL = 0.01;
%Background Parameters
k_b = 1.380*10^-23; %Boltzmann constant in Joules/Kelvin
n_b = 10^19; %background density in SI units (particles/m.^3)
T_b = 116.050; %background temperature in Kelvin
z_b = 1; %atomic number of background specie
m_b = 9.109*10^-31; %background mass in kg 
v_b = (2*k_b*T_b/m_b)^0.5; %background thermal velocity in m/s
eps_o = 8.85*10^-12; %permittivity of free space in Farad/m
n_o = 2*n_b; %initial number density in m^-3

%Target Specie Parameters
T_a = 1.1*T_b; %Target temperature in Kelvin
z_a = 1;
e = 1.602*10^-19; %charge in Coulombs
ln_Delt = 10; %Coulomb logarithm
m_a = 1.6726*10^-27; %target mass in kg
v_a = (2*k_b*T_a/m_a)^0.5; %target thermal velocity in m/s
L_ab = (e^2/(m_a*eps_o))^2; %Coefficient accounting for Coluomb force
nu_s = @(v) psi(v/v_b)*n_b*L_ab*(1 + m_a/m_b)./(2*pi.*v.^3); %Slowing down frequency in s^-1
nu_par = @(v) psi(v/v_b)*n_b*L_ab./(2*pi.*v.^3); %parallel diffusion frequency


%E = 1.0; %parallel Electric field

%% Setup the dimensions
% 
function ret = phi(x)
        ret = erf(x);
end

function ret = psi(x)     
        dphi_dx = 2./sqrt(pi) * exp(-x.^2);
        ret = 1./(2*x.^2) .* (phi(x) - x.*dphi_dx);   
        ix = find(abs(x)<1e-5); % catch singularity at boundary
        ret(ix) = 0;
 end
dim_v.domainMin = 0.1;
dim_v.domainMax = 10^4;
dim_v.init_cond_fn = @(x,p,t) n_o/(pi^3/2*v_a^3)*exp(-(x./v_a).^2);%.*(x == v_b);

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

% term V1 == 1/v^2 d/dv(v^3(m_a/(m_a + m_b))nu_s f))
% term V1 == g(v) q(v)      [mass, g(v) = 1/v^2,  BC N/A]
% q(v) == d/dv(g2(v)f(v))   [grad, g2(v) = v^3(m_a/(m_a + m_b))nu_s, BCL= N, BCR=N]

g1 = @(v,p,t,dat) 1./v.^2;
g2 = @(v,p,t,dat) v.^3*m_a.*nu_s(v)./(m_a + m_b);

pterm1 = MASS(g1);
pterm2  = GRAD(num_dims,g2,-1,'N','N');
termV_s = TERM_1D({pterm1,pterm2});
termV1   = TERM_ND(num_dims,{termV_s});

%%
% term V2 == 1/v^2 d/dv(v^4*0.5*nu_par*d/dv(f))
% term V2 == g(v) q(v)      [mass, g(v) = 1/v^2,  BC N/A]
% q(v) == d/dv(g2(v)r(v))   [grad, g2(v) = v^4*0.5*nu_par, BCL= D, BCR=D]
% r(v) = d/dv(g3(v)f)       [grad, g3(v) = 1, BCL=N, BCR=N]

g1 = @(v,p,t,dat) 1./v.^2;
g2 = @(v,p,t,dat) v.^4.*0.5.*nu_par(v);
g3 = @(v,p,t,dat) v.*0 + 1;

pterm1 = MASS(g1);
pterm2 = GRAD(num_dims,g2,-1,'D','D');
pterm3 = GRAD(num_dims,g3,+1,'N','N');
termV_par = TERM_1D({pterm1,pterm2,pterm3});
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
    @(x,p,t) (n_b/(pi^3/2*v_b^3)).*exp(-(x./v_b).^2), ...
    @(t,p) 1 
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