function pde = mirror_pitch
% PDE using the 1D Diffusion Equation for a case involving pitcha angle to axis \
% of the magnetic field line in a mirror configuration. This PDE is
% time dependent (although not all the terms are time dependent). The test
% particle is a hydrogen ion colliding with a background electron
% distribution that has a Maxwellian form. 
% PDE:
% 
% df/dt == nu_D/(2*sin(z)) d/dz ( sin(z) df/dz )
%
% nu_D is the deflection frequency
% Domain is [-pi,pi]
% Homogeneous Neumann boundary condition 
% Code will be added to equation involving velocity dimension
%
% Diffusion term is dealt with via LDG, i.e., splitting into two first
% order equations
%
% Run with
%
% explicit
% asgard(mirror_pitch);
%
% implicit
% asgard(mirror_pitch,'timestep_method','CN', 'deg', 3, 'lev', 3);
%

pde.CFL = 0.01;

%% Setup the dimensions
% 
% Here we setup a 1D problem 
    function ret = psi(x)
        
        phi = erf(x);
        dphi_dx = 2./sqrt(pi) * exp(-x.^2);

        ret = 1./(2*x.^2) .* (phi - x.*dphi_dx);
    end

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
nu_D = 10^9; %deflection frequency in s^-1 for v approximately equatl to v_b

%Initial parameters for target specie
n_o = 0.5*n_b; %initial number density for specie at specific velocity
v_o = 0.5*v_b; %initial known velocity

%function y = dd1(n)
% Our default value is 0
%y = 0; 

% The function is 1 only if the input is 0
%if n == 0
%    y = 1;
%end

%end
soln_z = @(z) n_b.*cos(z);
soln_t = @(t) t.*0 + 1;

%BCFunc = @(x) soln_x(x);
%BCFunc_t = @(t) soln_t(t);

% Domain is (a,b)

% The function is defined for the plane
% x = a and x = b
%BCL_fList = { ...
%    @(x,p,t) BCFunc(x), ... % replace x by a
%    @(t,p) BCFunc_t(t)
%    };

%BCR_fList = { ...
%    @(x,p,t) BCFunc(x), ... % replace x by b
%    @(t,p) BCFunc_t(t)
%    };

dim_z.name = 'z';
dim_z.domainMin = -pi;
dim_z.domainMax = pi;
dim_z.init_cond_fn = @(z,p,t) soln_z(z)*soln_t(t);


%%
% Add dimensions to the pde object
% Note that the order of the dimensions must be consistent with this across
% the remainder of this PDE.

pde.dimensions = {dim_z};
num_dims = numel(pde.dimensions);

%% Setup the terms of the PDE
%
% Here we have 1 term1, having only nDims=1 (x) operators.

%% 
% termC == nu_D/(2*sin(z))*d/dz sin(z)*df/dz
%
% becomes 
%
% termC == g1(z) q(z)        [mass, g1(p) = nu_D/(2*sin(z)),  BC N/A]
%   q(p) == d/dz g2(z) r(z)   [grad, g2(p) = sin(z), BCL=D,BCR=N]
%   r(p) == d/dp g3(z) f(z)   [grad, g3(p) = 1,      BCL=N,BCR=D]


g1 = @(z,p,t,dat) nu_D./(2*sin(z));
g2 = @(z,p,t,dat) sin(z);
g3 = @(z,p,t,dat) z.*0 + 1;
pterm1  = MASS(g1);
pterm2  = GRAD(num_dims,g2,+1,'N','N');
pterm3 = GRAD(num_dims,g3,-1,'D','D');
termC_z = TERM_1D({pterm1,pterm2,pterm3});
termC   = TERM_ND(num_dims,{termC_z});

%%
% Add terms to the pde object

pde.terms = {termC};

%% Construct some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

pde.params = params;

%%
% Sources

%s1x = @(x,p,t) -nu^2*cos(nu*x);
%s1t = @(t,p) exp(-2*nu^2*t);
%source1 = {s1x, s1t};

pde.sources = {};

%% Define the analytic solution (optional).
% This requires nDims+time function handles.

pde.analytic_solutions_1D = { ...
    @(z,p,t) soln_z(z), ...
    @(t,p) soln_t(t) 
    };

    function dt=set_dt(pde,CFL)
        
        dims = pde.dimensions;
        
        % for Diffusion equation: dt = C * dx^2
        
        lev = dims{1}.lev;
        dx = 1/2^lev;
        dt = CFL*dx^2;
        
    end

pde.set_dt = @set_dt;

end

%%
% Function to set time step