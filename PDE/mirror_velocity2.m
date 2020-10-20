function pde = mirror_velocity2(opts)
% Two-dimensional magnetic mirror from the FP paper - evolution of the ion velocity dependence
% of f in the presence of Coulomb collisions with background electrons
% 
% df/dt == nu_D/(2*sin(z)) d/dz ( sin(z) df/dz ) + 1/v^2 (d/dv(flux_v))
%a
% flux_v == v^3[(m_a/(m_a + m_b))nu_s f) + 0.5*nu_par*v*d/dv(f)]
%
% Run with
%
% asgard(@mirror_velocity2,'timestep_method','BE','case',3)

m_e = 9.109*10^-31; %electron mass in kg
m_H = 1.6726*10^-27; %hydrogen mass in kgs
m_D = 3.3443*10^-27; %deuterium mass in kgs
eps0 = 8.85*10^-12; %permittivity of free space in Farad/m
k_b = 1.380*10^-23; %Boltzmann constant in Joules/Kelvin
e = 1.602*10^-19; %charge in Coulombs
ln_delt = 12; %Coulomb logarithm

v_th = @(T_eV,m) sqrt(2*T_eV * e/m);

% species b: electrons in background
b.n = 4e19;
b.T_eV = 4;
b.Z = -1;
b.m = m_e;
b.vth = v_th(b.T_eV,b.m);

%species b2: deuterium in background
b2.n = 4e19;
b2.T_eV = 4;
b2.Z = 1;
b2.m = m_D;
b2.vth = v_th(b2.T_eV,b2.m);

% species a: electrons in beam
a.n = 4e19;
a.T_eV = 1e3;
a.Z = -1;
a.m = m_e;
a.vth = v_th(a.T_eV,a.m);

% L_ab = ln_delt*e^4/(m_a*eps0)^2; %Coefficient accounting for Coluomb force
switch opts.case_
    case 1
        T_eV_a = 0.05*T_eV_b; %Target temperature in Kelvin\
        offset = 10^6; %case with offset and change in Temperature
    case 2
        T_eV_a = 0.05*T_eV_b;
        offset = 0; %case with no offset but change in Temperature
    case 3
        a.T_eV = 1e3;
        offset = 10^7; %case with offset and no change in Temperature
end

% T_a = T_eV_a*11606; %converting to Kelvin
% T_b = T_eV_b*11606; %converting to Kelvin
% nu_s = @(v) psi(v./v_th(T_b,m_b)).*n_b*L_ab*(1 + m_a/m_b)./(2*pi*v_th(T_b,m_b).^3.*v./v_th(T_b,m_b)); %Slowing down frequency in s^-1; 


x = @(v,vth) v./vth; 
nu_ab0 = @(a,b) b.n * e^4 * a.Z^2 * b.Z^2 * ln_delt / (2*pi*eps0^2*a.m^2*b.vth^3); %scaling coefficient
nu_s = @(v,a,b) nu_ab0(a,b) .* (1+a.m/b.m) .* psi(x(v,b.vth)) ./ x(v,b.vth); %slowing down frequency
nu_par = @(v,a,b) nu_ab0(a,b).*(psi(x(v,b.vth))./(x(v,b.vth).^3)); %parallel diffusion frequency
nu_D =  @(v,a,b) nu_ab0(a,b).*(phi_f(x(v,b.vth)) - psi(x(v,b.vth)))./(x(v,b.vth).^3); %deflection frequency in s^-1

maxwell = @(v,x,y) a.n/(pi^3/2.*y^3).*exp(-((v-x)/y).^2);
gauss = @(v,x) a.n/(sqrt(2*pi)*x)*exp(-0.5*((v - x)/x).^2);
% norm = 3.749;
init_func = @(v) maxwell(v,0,v_th(b.T_eV,a.m)) + maxwell(v,5*10^6, 10^6);
pitch_z = @(z) z.*0 + 1;
pitch_t = @(t) exp(-nu_D(v_th(b.T_eV,a.m),a,b).*t);

v_ = 10.^[1:.1:7];
%loglog(0.5*a.m*v_.^2/e,nu_D(v_,a,a))
hold on
%loglog(0.5*a.m*v_.^2/e,nu_D(v_,a,b))
xlim([0.1,1e2]);
ylim([1e4,1e11]);

BCFunc = @(v) init_func(v);


% Domain is (a,b)

% The function is defined for the plane
% x = a and x = b
BCL_fList = { ...
    @(v,p,t) v.*0, ... 
    @(z,p,t) pitch_z(z), ...
    @(t,p) pitch_t(t)
    };

BCR_fList = { ...
    @(v,p,t) BCFunc(v), ... % 
    @(z,p,t) pitch_z(z), ...
    @(t,p) pitch_t(t)
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

%% Define the dimensions

dim_v = DIMENSION(0,3e7);
dim_v.name = 'v';
dim_v.init_cond_fn = @(v,p,t) init_func(v);
dim_v.jacobian = @(v,p,t) 2.*pi.*v.^2;

dim_z = DIMENSION(0,pi);
dim_z.name = 'z';
dim_z.init_cond_fn = @(z,p,t) pitch_z(z)*pitch_t(t);
dim_z.jacobian = @(z,p,t) sin(z);

dimensions = {dim_v, dim_z};
num_dims = numel(dimensions);

%% Define the terms of the PDE

%% 
% termC == nu_D/(2*sin(z))*d/dz sin(z)*df/dz
%
% becomes 
%
% termC == g1(v) g2(z) q(z)   [mass, g1(p) = nu_D(v), g2(z) = 1/(2sin(z))  BC N/A]
%   q(z) == d/dz g3(z) r(z)   [grad, g3(z) =  sin(z), BCL=D,BCR=D]
%   r(z) == d/dp g4(z) f(z)   [grad, g3(p) = 1,      BCL=N,BCR=N]


g1 = @(v,p,t,dat) nu_D(v,a,b);
g2 = @(z,p,t,dat) 1./(2*sin(z));
g3 = @(z,p,t,dat) sin(z);
g4 = @(z,p,t,dat) z.*0 + 1;
pterm1  = MASS(g1);
pterm2  = MASS(g2);
pterm3  = GRAD(num_dims,g3,+1,'D','D');
pterm4 = GRAD(num_dims,g4,-1,'N', 'N');
termC_z = TERM_1D({pterm1,pterm2,pterm3,pterm4});
termC   = TERM_ND(num_dims,{termC_z,[]});

% term V1 == 1/v^2 d/dv(v^3(m_a/(m_a + m_b))nu_s f))
% term V1 == g(v) q(v)      [mass, g(v) = 1/v^2,  BC N/A]
% q(v) == d/dv(g2(v)f(v))   [grad, g2(v) = v^3(m_a/(m_a + m_b))nu_s, BCL= N, BCR=N]

g1 = @(v,p,t,dat) 1./v.^2;
g2 = @(v,p,t,dat) v.^3*a.m.*nu_s(v,a,b)./(a.m + b.m);

pterm1 = MASS(g1);
pterm2  = GRAD(num_dims,g2,-1,'N','D', BCL_fList, BCR_fList);
termV_s = TERM_1D({pterm1,pterm2});
termV1   = TERM_ND(num_dims,{termV_s,[]});

%%
% term V2 == 1/v^2 d/dv(v^4*0.5*nu_par*d/dv(f))
% term V2 == g(v) q(v)      [mass, g(v) = 1/v^2,  BC N/A]
% q(v) == d/dv(g2(v)r(v))   [grad, g2(v) = v^4*0.5*nu_par, BCL= D, BCR=D]
% r(v) = d/dv(g3(v)f)       [grad, g3(v) = 1, BCL=N, BCR=N]

g1 = @(v,p,t,dat) 1./v.^2;
g2 = @(v,p,t,dat) v.^4.*0.5.*nu_par(v,a,b);
g3 = @(v,p,t,dat) v.*0 + 1;

pterm1 = MASS(g1);
pterm2 = GRAD(num_dims,g2,-1,'D','N');
pterm3 = GRAD(num_dims,g3,+1,'N','D', BCL_fList, BCR_fList);
termV_par = TERM_1D({pterm1,pterm2,pterm3});
termV2   = TERM_ND(num_dims,{termV_par,[]});

terms = {termV1,termV2, termC};

%% Define some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

%% Define sources

sources = {};

%% Define the analytic solution (optional).
% This requires nDims+time function handles.

    function ans = soln(v,p,t)
        ans = a.n/(pi^3/2.*v_th(b.T_eV,a.m).^3).*exp(-(v./v_th(b.T_eV,a.m)).^2);        
        if isfield(p,'norm_fac')
            ans = ans .* p.norm_fac;
        end
    end
    soln_handle = @soln;
    
analytic_solutions_1D = { ...    
    @(v,p,t) soln_handle(v,p,t), ...
    @(z,p,t) pitch_z(z), ...
    @(t,p) t.*0 + 1; %pitch_t(t)
    };

%% Define function to set time step
    function dt=set_dt(pde,CFL)    
        Lmax = pde.dimensions{1}.max;
        LevX = pde.dimensions{1}.lev;
        dt = Lmax/2^LevX*CFL;
    end

%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,analytic_solutions_1D);

end