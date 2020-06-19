function pde = diffusion1_mirror
% PDE using the 1D Diffusion Equation for a case involving velocity
% parallel to the magnetic field line in a mirror configuration. This PDE is
% time dependent (although not all the terms are time dependent). The test
% particle is a hydrogen ion colliding with a background electron
% distribution that has a Maxwellian form. 
% PDE:
% 
% df/dt = d/dv_|| (-D(v_||) df_a /dv_||
% Domain is [-1,1]
% Homogeneous Neumann boundary condition 
% Code will be added to equation involving perpendicular velocity term 
%
% Diffusion terms are dealt with via LDG, i.e., splitting into two first
% order equations
%
% Run with
%
% explicit
% asgard(diffusion1);
%
% implicit
% asgard(diffusion1_mirror,'timestep_method','CN');
%
%Louis: With small nu, stability is maintained for longer times
% asgard(diffusion1,'lev',3,'deg',4,'timestep_method','BE', 'dt',0.05,'num_steps',20)
pde.CFL = 0.01;

%% Setup the dimensions
% 
% Here we setup a 1D problem 
    function ret = psi(x,t)
        
        phi = erf(x);
        dphi_dx = 2./sqrt(pi) * exp(-x.^2);

        ret = 1./(2*x.^2) .* (phi - x.*dphi_dx);
    end

%Background Parameters
n_b = 10^15; %background density in SI units (particles/m.^3)
T_b = 1000; %background temperature in eV
z_b = 1; %atomic number of background specie
m_b = 0.0005485; %background mass in amu 
v_b = sqrt(2*T_b/m_b)*9822.694; %background velocity in m/s

%Target Specie Parameters
z_a = 1;
e = 1.602*10^-19; %charge in Coulombs
ln_Delt = 10; %Coulomb logarithm
m_a = 1.6726*10^-19; %target mass in kg
Gamma_a = ln_Delt*e^4*z_a^2/m_a^2;

%Initial parameters for target specie
n_o = 0.5*n_b; %initial number density for specie at specific velocity
v_o = 0.5*v_b; %initial known velocity

function y = dd1(n)
% Our default value is 0
y = 0; 

% The function is 1 only if the input is 0
if n == 0
    y = 1;
end

end
soln_v = @(x) x.*0 + n_o.*exp(-(x-v_o).^2/(0.1*v_o).^2);
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

dim_x.name = 'v_par';
dim_x.domainMin = 0.001;
dim_x.domainMax = v_b;
dim_x.init_cond_fn = @(x,p,t) soln_v(x)*soln_t(t);


%%
% Add dimensions to the pde object
% Note that the order of the dimensions must be consistent with this across
% the remainder of this PDE.

pde.dimensions = {dim_x};
num_dims = numel(pde.dimensions);

%% Setup the terms of the PDE
%
% Here we have 1 term, with each term having nDims (x and y) operators.

%% 
% Setup the d^2_dx^2 term

% term1
%
% eq1 :  df/dt   == d/dx g1(x) q(x,y)   [grad,g1(x)=1, BCL=N, BCR=D]
% eq2 :   q(x,y) == d/dx g2(x) f(x,y)   [grad,g2(x)=1, BCL=D, BCR=D]
% coeff_mat = mat1 * mat2

g1 = @(x,p,t,dat) x.*0 + 1;
% g2 = @(x,p,t,dat) -Gamma_a*n_b*(z_b)^2.*psi(abs(x)/v_b,t)./abs(x);
g2 = @(x,p,t,dat) +1e14;

pterm1 = GRAD(num_dims,g1,+1,'D','D');
pterm2 = GRAD(num_dims,g2,-1,'N','N');

term1_x = TERM_1D({pterm1,pterm2});
term1   = TERM_ND(num_dims,{term1_x});

%%
% Add terms to the pde object

 pde.terms = {term1};

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
    @(x,p,t) soln_v(x), ...
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