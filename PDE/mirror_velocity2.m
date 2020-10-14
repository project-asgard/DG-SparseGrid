function pde = mirror_velocity2
% Two-dimensional magnetic mirror from the FP paper - evolution of the ion velocity dependence
% of f in the presence of Coulomb collisions with background electrons
% 
% df/dt == nu_D/(2*sin(z)) d/dz ( sin(z) df/dz ) + 1/v^2 (d/dv(flux_v))
%a
% flux_v == v^3[(m_a/(m_a + m_b))nu_s f) + 0.5*nu_par*v*d/dv(f)]
%
% Run with
%
% asgard(mirror_velocity2,'timestep_method','BE')

params = mirror_parameters();
pde.params = params;

test = 'c';
pde.CFL = 0.01;

switch test
    case 'a'
        params.a.T_eV = 0.05*params.b.T_eV; %Target temperature in Kelvin\
        offset = 10^6; %case with offset and change in Temperature
    case 'b'
        params.a.T_eV = 0.05*params.b.T_eV;
        offset = 0; %case with no offset but change in Temperature
    case 'c'
        params.a.T_eV = 1e3;
        offset = 10^7; %case with offset and no change in Temperature
end

BCFunc = @(v,p,t) p.init_cond_v(v);

% Domain is (a,b)

% The function is defined for the plane
% x = a and x = b
BCL_fList = { ...
    @(v,p,t) v.*0, ... 
    @(z,p,t) p.init_cond_z(z), ...
    @(t,p) p.init_cond_t(t)
    };

BCR_fList = { ...
    @(v,p,t) BCFunc(v,p,t), ... % 
    @(z,p,t) p.init_cond_z(z), ...
    @(t,p) p.init_cond_t(t)
    };

%% Setup the dimensions
% 
dim_v.name = 'v';
dim_v.domainMin = 0;
dim_v.domainMax = 3*10^7;
dim_v.init_cond_fn = @(v,p,t) p.init_cond_v(v);
dim_v.jacobian = @(v,p,t) 2.*pi.*v.^2;

dim_z.name = 'z';
dim_z.domainMin = 0;
dim_z.domainMax = pi;
dim_z.init_cond_fn = @(z,p,t) p.init_cond_z(z)*p.init_cond_t(t);
dim_z.jacobian = @(z,p,t) sin(z);

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
% termC == g1(v) g2(z) q(z)   [mass, g1(p) = nu_D(v), g2(z) = 1/(2sin(z))  BC N/A]
%   q(z) == d/dz g3(z) r(z)   [grad, g3(z) =  sin(z), BCL=D,BCR=D]
%   r(z) == d/dp g4(z) f(z)   [grad, g3(p) = 1,      BCL=N,BCR=N]

g1 = @(v,p,t,dat) p.nu_D(v,p.a,p.b);
g2 = @(z,p,t,dat) 1./(2*sin(z));
g3 = @(z,p,t,dat) sin(z);
g4 = @(z,p,t,dat) z.*0 + 1;
pterm1  = MASS(g1);
pterm2  = MASS(g2);
pterm3  = GRAD(num_dims,g3,+1,'D','D');
pterm4  = GRAD(num_dims,g4,-1,'N', 'N');
termC_z = TERM_1D({pterm1,pterm2,pterm3,pterm4});
termC   = TERM_ND(num_dims,{termC_z,[]});

% term V1 == 1/v^2 d/dv(v^3(m_a/(m_a + m_b))nu_s f))
% term V1 == g(v) q(v)      [mass, g(v) = 1/v^2,  BC N/A]
% q(v) == d/dv(g2(v)f(v))   [grad, g2(v) = v^3(m_a/(m_a + m_b))nu_s, BCL= N, BCR=N]

g1 = @(v,p,t,dat) 1./v.^2;
g2 = @(v,p,t,dat) v.^3*p.a.m.*p.nu_s(v,p.a,p.b)./(p.a.m + p.b.m);

pterm1  = MASS(g1);
pterm2  = GRAD(num_dims,g2,-1,'N','D', BCL_fList, BCR_fList);
termV_s = TERM_1D({pterm1,pterm2});
termV1  = TERM_ND(num_dims,{termV_s,[]});

%%
% term V2 == 1/v^2 d/dv(v^4*0.5*nu_par*d/dv(f))
% term V2 == g(v) q(v)      [mass, g(v) = 1/v^2,  BC N/A]
% q(v) == d/dv(g2(v)r(v))   [grad, g2(v) = v^4*0.5*nu_par, BCL= D, BCR=D]
% r(v) = d/dv(g3(v)f)       [grad, g3(v) = 1, BCL=N, BCR=N]

g1 = @(v,p,t,dat) 1./v.^2;
g2 = @(v,p,t,dat) v.^4.*0.5.*p.nu_par(v,p.a,p.b);
g3 = @(v,p,t,dat) v.*0 + 1;

pterm1      = MASS(g1);
pterm2      = GRAD(num_dims,g2,-1,'D','N');
pterm3      = GRAD(num_dims,g3,+1,'N','D', BCL_fList, BCR_fList);
termV_par   = TERM_1D({pterm1,pterm2,pterm3});
termV2      = TERM_ND(num_dims,{termV_par,[]});

%%
% Add terms to the pde object

pde.terms = {termV1,termV2, termC};

%% Add an arbitrary number of sources to the RHS of the PDE
% Each source term must have nDims + 1

%%
% Add sources to the pde data structure
pde.sources = {};

%% Define the analytic solution (optional).
% This requires nDims+time function handles.

pde.analytic_solutions_1D = { ...    
    @(v,p,t) p.analytic_solution_v(v,p,t), ...
    @(z,p,t) p.init_cond_z(z), ...
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