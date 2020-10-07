function pde = mirror_velocity
% One-dimensional magnetic mirror from the FP paper - evolution of the ion velocity dependence
% of f in the presence of Coulomb collisions with background electrons
% 
% df/dt == 1/v^2 (d/dv(flux_v))
%a
% flux_v == v^3[(m_a/(m_a + m_b))nu_s f) + 0.5*nu_par*v*d/dv(f)]
%
% Run with
%
% asgard(mirror_velocity,'timestep_method','BE')

test = 'c';

pde.CFL = 0.01;

m_e = 9.109*10^-31; %electron mass in kg
m_D = 3.3443*10^-27; %deuterium mass in kgs
m_H = 1.6726*10^-27; %hydrogen mass in kgs
eps_o = 8.85*10^-12; %permittivity of free space in Farad/m
k_b = 1.380*10^-23; %Boltzmann constant in Joules/Kelvin
e = 1.602*10^-19; %charge in Coulombs
ln_delt = 10; %Coulomb logarithm

%Background Species parameter
n_b = 4*10^19; %background density in SI units (particles/m.^3)
T_eV_b = 4; %background temperature in eV
z_b = 1; %atomic number of background
m_b = m_D; %background mass

%Target Specie Parameters
n_a = n_b;
z_a = 1;
m_a = m_e;%target species

v_th = @(T,m) (2*k_b*T./m).^0.5; %thermal velocity function
L_ab = ln_delt*e^4/(m_a*eps_o)^2; %Coefficient accounting for Coluomb force

switch test
    case 'a'
        T_eV_a = 0.05*T_eV_b; %Target temperature in Kelvin\
        offset = 10^6; %case with offset and change in Temperature
    case 'b'
        T_eV_a = 0.05*T_eV_b;
        offset = 0; %case with no offset but change in Temperature
    case 'c'
        T_eV_a = 1000;
        T_a = T_eV_a*11606; %converting to Kelvin
        offset = v_th(T_a,m_a); %case with offset and no change in Temperature
end 
T_b = T_eV_b*11606; %converting to Kelvin
nu_ab0 = n_b*e^4*z_a^2*z_b^2*ln_delt/(2*pi*eps_o^2*m_a^2*v_th(T_b,m_b)^3);
nu_s = @(v) nu_ab0.*(1 + m_a/m_b).*psi(v./v_th(T_b,m_b))./(v./v_th(T_b,m_b)); %Slowing down frequency in s^-1; %parallel diffusion frequency
loglog(0.5*m_a*logspace(-1,7).^2/e, nu_s(logspace(-1,7)))
%ylim([10^4, 10^11]);
xlim([10^-1, 10^2]);
nu_par = @(v) psi(v./v_th(T_b,m_b)).*n_b*L_ab./(2*pi.*v.^3);
maxwell = @(v,x,y) n_a/(pi^3/2.*y^3).*exp(-((v-x)/y).^2);
gauss = @(v,x) n_a/(sqrt(2*pi)*x)*exp(-0.5*((v - x)/x).^2);
norm = v_th(T_b,m_a)*(sqrt(pi)*(v_th(T_a,m_a)^2 + 2*offset^2)*phi((offset - 0.01)/v_th(T_a,m_a)) + 2*v_th(T_a,m_a)*(0.01 + offset)*exp(-(0.01 - offset)^2/v_th(T_a,m_a)^2) + sqrt(pi)*(v_th(T_a,m_a)^2 +2*offset^2)*phi((10^7 - offset)/v_th(T_a,m_a)) - 2*v_th(T_a,m_a)*(10^7 + offset)*exp(-(10^7 - offset)^2/v_th(T_a,m_a)^2))/(v_th(T_a,m_a)^2*(2*0.01*exp(-0.01^2/v_th(T_b,m_a)^2) - 2*10^7*exp(-10^14/v_th(T_b,m_a)^2) + sqrt(pi)*v_th(T_b,m_a)*(phi(10^7/v_th(T_b,m_a)) - phi(0.01/v_th(T_b,m_a)))));
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
%soln_t = @(t) exp(-m_a*nu_s.*t/(m_a + m_b));

dim_v.domainMin = 0;
dim_v.domainMax = 3*10^7;
dim_v.jacobian = @(v,p,t) 4.*pi.*v.^2;

%steady state analytic solution
soln_v = @(v) n_a*(2*(m_a/(m_a + m_b))*nu_s./nu_par(v_th(T_b,m_a))).*(1/(4*pi)).*(1./v.^3);
%normalized coordinates by thermal velocity
%v_n_max = dim_v.domainMax/v_th(T_b,m_b);
%v_n_min = dim_v.domainMin/v_th(T_b,m_b);
%initial condition
var = 10^6;
dim_v.init_cond_fn = @(v,p,t) maxwell(v,v_th(T_a,m_a), var);%.*(x == v_b);

BCFunc = @(v) maxwell(v,v_th(T_b,m_a), v_th(T_b,m_a));

% Domain is (a,b)

% The function is defined for the plane
% x = a and x = b
BCL_fList = { ...
    @(v,p,t) v.*0, ... 
    @(t,p) t.*0 + 1
    };

BCR_fList = { ...
    @(v,p,t) BCFunc(v), ... % 
    @(t,p) t.*0 + 1
    };

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
g2 = @(v,p,t,dat) v.^3.*nu_s(v).*m_a/(m_a + m_b);

pterm1 = MASS(g1);
pterm2  = GRAD(num_dims,g2,-1,'N','D', BCL_fList, BCR_fList);
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
pterm2 = GRAD(num_dims,g2,+1,'D','N');
pterm3 = GRAD(num_dims,g3,-1,'N','D', BCL_fList, BCR_fList);
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

source = { ...
    @(v,p,t) (m_a.*nu_s/(m_a+m_b))*(1-m_a*nu_s/(m_a + m_b)*(4*v+3*t)),   ...   % s1v
    @(t,p) t.*0 + 1 ...   % s1t
    };

%%
% Add sources to the pde data structure
pde.sources = {};

%% Define the analytic solution (optional).
% This requires nDims+time function handles.
    function ans = soln(v,p,t)
        ans = n_a/(pi^3/2.*v_th(T_b,m_a).^3).*exp(-(v./v_th(T_b,m_a)).^2);        
        if isfield(p,'norm_fac')
            ans = ans .* p.norm_fac;
        end
    end
    soln_handle = @soln;


pde.analytic_solutions_1D = { ...    
    @(v,p,t) soln_handle(v,p,t); ...
    @(t,p) t.*0 + 1; 
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