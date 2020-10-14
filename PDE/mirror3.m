function pde = mirror3
% Three-dimensional magnetic mirror from the FP paper - evolution of the ion velocity and spatial dependence
% of f in the presence of Coulomb collisions with background electrons
% 
% df/dt == -dB/ds f - vcos(z)df/dz + nu_D/(2*sin(z)) d/dz ( sin(z) df/dz ) + 1/v^2 (d/dv(flux_v))
%a
% flux_v == v^3[(m_a/(m_a + m_b))nu_s f) + 0.5*nu_par*v*d/dv(f)]
%
% Run with
%
% asgard(mirror3,'timestep_method','BE')

test = 'c';
pde.CFL = 0.01;

 m_e = []; %electron mass in kg
 m_H = []; %hydrogen mass in kgs
 m_D = []; %deuterium mass in kgs
 eps0 = []; %permittivity of free space in Farad/m
 k_b = []; %Boltzmann constant in Joules/Kelvin
 e = []; %charge in Coulombs
 ln_delt = []; %Coulomb logarithm

 v_th = [];

% species b: electrons in background
 b.n = [];
% 
% %species b2: deuterium in background
 b2.n = [];

% % species a: electrons in beam
 a.n = [];
 
switch test
    case 'a'
        a.T_eV = 0.05*T_eV_b; %Target temperature in Kelvin\
        offset = 10^6; %case with offset and change in Temperature
    case 'b'
        a.T_eV = 0.05*T_eV_b;
        offset = 0; %case with no offset but change in Temperature
    case 'c'
        a.T_eV = 1e3;
        offset = 10^7; %case with offset and no change in Temperature
end

 x = []; 
 nu_ab0 = []; %scaling coefficient
 nu_s = []; %slowing down frequency
 nu_par = []; %parallel diffusion frequency
 nu_D = [];
 
maxwell = [];
gauss = [];

mirror_common

% norm = 3.749;
init_func_v = @(v) maxwell(v,v_th(a.T_eV,a.m), 10^6);
init_func_z = @(z) z.*0 + 1;

v_ = 10.^[-1:.1:7];
%loglog(0.5*a.m*v_.^2/e,nu_D(v_,a,b))
hold on
%loglog(0.5*a.m*v_.^2/e,nu_s(v_,a,b))
%loglog(0.5*a.m*v_.^2/e,nu_par(v_,a,b))
xlim([0.1,1e2]);
ylim([1e4,1e11]);

BCFunc = @(v) init_func(v);

B_func = @(s) sin(s); %magnetic field as a function of spatial coordinate
dB_ds = @(s) cos(s); %derivative of magnetic field
v_test = 500; %test value for velocity in spatial advection coefficient
z_test = pi/2 - 1e-6; %test value for pitch angle in spatial advection coefficient
A = v_test*cos(z_test); %spatial advection coefficient
%Expected solution for velocity space, with norm to account for particle
%conservation
%Expected solution in pitch angle
pitch_z = @(z) n_a.*cos(z);
pitch_t = @(t) exp(-nu_D*t);
%Expected solution for spatial dimension
space_s = @(s) exp(s);
space_t = @(t) exp(-2*A*t);
%Boundary for velocity space
BCFunc = @(v) maxwell_func(v,T_b,m_a);


% Domain is (a,b)

% The function is defined for the plane
% x = a and x = b
BCL_fList = { ...
    @(v,p,t) v.*0, ... 
    @(z,p,t) pitch_z(z), ...
    @(s,p,t) space_s(s), ...
    @(t,p) pitch_t(t).*space_t(t)
    };

BCR_fList = { ...
    @(v,p,t) BCFunc(v), ... % 
    @(z,p,t) pitch_z(z), ...
    @(s,p,t) space_s(s), ...
    @(t,p) pitch_t(t).*space_t(t)
    };

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
dim_v.init_cond_fn = @(v,p,t) maxwell_func(v,T_a,m_a);

dim_z.name = 'z';
dim_z.domainMin = -pi;
dim_z.domainMax = pi;
dim_z.init_cond_fn = @(z,p,t) pitch_z(z)*pitch_t(t);

dim_s.domainMin = 0;
dim_s.domainMax = 5;
dim_s.init_cond_fn = @(s,p,t) space_s(s);

%%
% Add dimensions to the pde object
% Note that the order of the dimensions must be consistent with this across
% the remainder of this PDE.

pde.dimensions = {dim_v, dim_z, dim_s};
num_dims = numel(pde.dimensions);

%% Setup the terms of the PDE

%% 
%% Advection  term
% termS1 == -vcos(z)*df/ds
% termS1 == q(v)*r(z)*w(s)
% q(v) == g1(v)  [mass, g1(p) = v,  BC N/A]
% r(z) == g2(z) [mass, g2(z) = cos(z),  BC N/A]
% w(s) == d/ds g3(s) f [grad, g3(s) = -1, BCL= D, BCR=D]
g1 = @(v,p,t,dat) v;
g2 = @(z,p,t,dat) cos(z);
g3 = @(s,p,t,dat) s.*0 - 1;

pterm1 = MASS(g1);
pterm2 = MASS(g2);
pterm3 = GRAD(num_dims,g3,-1,'D','D', BCL_fList, BCR_fList);

term1_s = TERM_1D({pterm1,pterm2,pterm3});
termS1   = TERM_ND(num_dims,{term1_s,[],[]});

%%Mass term
%termS2 == -vcos(z)dB/ds f
% termS1 == q(v)*r(z)*w(s)
% q(v) == g1(v)  [mass, g1(p) = v,  BC N/A]
% r(z) == g2(z) [mass, g2(z) = cos(z),  BC N/A]
% w(s) == g3(s) f [mass, g3(s) = -dB/ds/B, BCL= D, BCR=D]

g1 = @(v,p,t,dat) v;
g2 = @(z,p,t,dat) cos(z);
g3 = @(s,p,t,dat) -dB_ds(s)./B_func(s);

pterm1 = MASS(g1);
pterm2 = MASS(g2);
pterm3 = MASS(g3);

termB1 = TERM_1D({pterm1,pterm2,pterm3});

termS2 = TERM_ND(num_dims,{termB1,[],[]});

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
termC   = TERM_ND(num_dims,{termC_z,[],[]});

% term V1 == 1/v^2 d/dv(v^3(m_a/(m_a + m_b))nu_s f))
% term V1 == g(v) q(v)      [mass, g(v) = 1/v^2,  BC N/A]
% q(v) == d/dv(g2(v)f(v))   [grad, g2(v) = v^3(m_a/(m_a + m_b))nu_s, BCL= N, BCR=D]

g1 = @(v,p,t,dat) 1./v.^2;
g2 = @(v,p,t,dat) v.^3*m_a*nu_s/(m_a + m_b);

pterm1 = MASS(g1);
pterm2  = GRAD(num_dims,g2,-1,'N','D', BCL_fList, BCR_fList);
termV_s = TERM_1D({pterm1,pterm2});
termV1   = TERM_ND(num_dims,{termV_s,[],[]});

%%
% term V2 == 1/v^2 d/dv(v^4*0.5*nu_par*d/dv(f))
% term V2 == g(v) q(v)      [mass, g(v) = 1/v^2,  BC N/A]
% q(v) == d/dv(g2(v)r(v))   [grad, g2(v) = v^4*0.5*nu_par, BCL= D, BCR=D]
% r(v) = d/dv(g3(v)f)       [grad, g3(v) = 1, BCL=N, BCR=N]

g1 = @(v,p,t,dat) 1./v.^2;
g2 = @(v,p,t,dat) v.^4.*0.5.*nu_par(v);
g3 = @(v,p,t,dat) v.*0 + 1;

pterm1 = MASS(g1);
pterm2 = GRAD(num_dims,g2,-1,'D','N');
pterm3 = GRAD(num_dims,g3,+1,'N','D', BCL_fList, BCR_fList);
termV_par = TERM_1D({pterm1,pterm2,pterm3});
termV2   = TERM_ND(num_dims,{termV_par,[],[]});

%%
% Add terms to the pde object

pde.terms = {termV1,termV2, termC,termS1,termS2};

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
    @(v,p,t) norm*n_a/(pi^3/2.*v_th(T_b,m_a).^3).*exp(-(v./v_th(T_b,m_a)).^2), ...
    @(z,p,t) pitch_z(z), ...
    @(s,p,t) space_s(s), ...
    @(t,p) t.*0 + 1; %pitch_t(t)
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