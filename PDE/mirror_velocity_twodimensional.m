function pde = mirror_velocity_twodimensional
% Two-dimensional magnetic mirror from the FP paper - evolution of the ion velocity dependence
% of f in the presence of Coulomb collisions with background electrons
% 
% df/dt == nu_D/(2*sin(z)) d/dz ( sin(z) df/dz ) + 1/v^2 (d/dv(flux_v))
%a
% flux_v == v^3[(m_a/(m_a + m_b))nu_s f) + 0.5*nu_par*v*d/dv(f)]
%
% Run with
%
% asgard(mirror_velocity_twodimensional,'timestep_method','BE')

pde.CFL = 0.01;

m_e = 9.109*10^-31; %electron mass in kg
m_H = 1.6726*10^-27; %hydrogen mass in kgs
eps_o = 8.85*10^-12; %permittivity of free space in Farad/m
k_b = 1.380*10^-23; %Boltzmann constant in Joules/Kelvin
e = 1.602*10^-19; %charge in Coulombs
ln_delt = 10; %Coulomb logarithm

%Background Species parameter
n_b = 10^19; %background density in SI units (particles/m.^3)
T_eV_b = 10; %background temperature in eV
z_b = 1; %atomic number of background
m_b = m_e; %background mass

%Target Specie Parameters
n_a = n_b;
T_eV_a = T_eV_b; %Target temperature in Kelvin
z_a = 1;
m_a = m_H;%target species

T_a = T_eV_a*11606; %converting to Kelvin
T_b = T_eV_b*11606; %converting to Kelvin

v_th = @(T,m) (2*k_b*T./m).^0.5; %thermal velocity function
L_ab = ln_delt*e^4/(m_a*eps_o)^2; %Coefficient accounting for Coluomb force
nu_s = @(v) psi(v./v_th(T_b,m_b)).*n_b*L_ab*(1 + m_a/m_b)./(2*pi*v_th(T_b,m_b).^3.*v./v_th(T_b,m_b)); %Slowing down frequency in s^-1
nu_par = @(v) psi(v./v_th(T_b,m_b)).*n_b*L_ab./(2*pi.*v.^3); %parallel diffusion frequency
nu_D = 10^4; % @(v) n_b*L_ab.*(phi_f(v./v_th(T_b,m_b)) - psi(v./v_th(T_b,m_b)))./(4*pi.*v.^3); %deflection frequency in s^-1
offset = 10^6;
maxwell_func = @(v,T,m) n_b/(pi^3/2.*v_th(T,m).^3).*exp(-(v./v_th(T,m)).^2);
pitch_z = @(z) n_a.*cos(z);
pitch_t = @(t) exp(-nu_D*t);

%E = 1.0; %parallel Electric field

function ret = phi(x)
        ret = erf(x);
end
function ret = phi_f(x)
        ret = (x + 1./(2*x)).*erf(x) + exp(-x.^2)./sqrt(pi);     
end
function ret = psi(x)     
        dphi_dx = 2./sqrt(pi) * exp(-x.^2);
        ret = 1./(2*x.^2) .* (phi(x) - x.*dphi_dx);   
        ix = find(abs(x)<1e-5); % catch singularity at boundary
        ret(ix) = 0;
end

%% Setup the dimensions
% 
dim_v.name = 'v';
dim_v.domainMin = 0.01;
dim_v.domainMax = 10^7;
dim_v.init_cond_fn = @(v,p,t) n_b/(pi^3/2.*v_th(T_b,m_b).^3).*(v == v_th(T_b,m_a));

dim_z.name = 'z';
dim_z.domainMin = -pi;
dim_z.domainMax = pi;
dim_z.init_cond_fn = @(z,p,t) pitch_z(z)*pitch_t(t);

%%
% Add dimensions to the pde object
% Note that the order of the dimensions must be consistent with this across
% the remainder of this PDE.

pde.dimensions = {dim_v, dim_z};
num_dims = numel(pde.dimensions);

%% Setup the terms of the PDE

%% 
% termC == nu_D/(2*sin(z))*d/dz sin(z)*df/dz
%
% becomes 
%
% termC == g1(z) q(z)        [mass, g1(p) = nu_D/(2*sin(z)),  BC N/A]
%   q(p) == d/dz g2(z) r(z)   [grad, g2(p) = sin(z), BCL=N,BCR=D]
%   r(p) == d/dp g3(z) f(z)   [grad, g3(p) = 1,      BCL=D,BCR=N]


g1 = @(z,p,t,dat) nu_D./(2*sin(z));
g2 = @(z,p,t,dat) sin(z);
g3 = @(z,p,t,dat) z.*0 + 1;
pterm1  = MASS(g1);
pterm2  = GRAD(num_dims,g2,+1,'D','D');
pterm3 = GRAD(num_dims,g3,0,'N', 'N');
termC_z = TERM_1D({pterm1,pterm2,pterm3});
termC   = TERM_ND(num_dims,{termC_z,[]});

% term V1 == 1/v^2 d/dv(v^3(m_a/(m_a + m_b))nu_s f))
% term V1 == g(v) q(v)      [mass, g(v) = 1/v^2,  BC N/A]
% q(v) == d/dv(g2(v)f(v))   [grad, g2(v) = v^3(m_a/(m_a + m_b))nu_s, BCL= N, BCR=N]

g1 = @(v,p,t,dat) 1./v.^2;
g2 = @(v,p,t,dat) v.^3*m_a.*nu_s(v)./(m_a + m_b);

pterm1 = MASS(g1);
pterm2  = GRAD(num_dims,g2,-1,'N','N');
termV_s = TERM_1D({pterm1,pterm2});
termV1   = TERM_ND(num_dims,{termV_s,[]});

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
termV2   = TERM_ND(num_dims,{termV_par,[]});

%%
% Add terms to the pde object

pde.terms = {termV1,termV2, termC};

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
    @(v,p,t) maxwell_func(v, T_b, m_a), ...
    @(z,p,t) pitch_z(z)*pitch_t(t), ...
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