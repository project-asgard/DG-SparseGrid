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
% asgard(@mirror_velocity2,'timestep_method','BE','case',3,'dt',1e-6)

params = mirror_parameters();

switch opts.case_
    case 1 
        params.a.T_eV = 0.05*params.b.T_eV; %Target temperature in Kelvin\
        offset = 10^6; %case with offset and change in Temperature
    case 2 
        params.a.T_eV = 0.05*params.b.T_eV;
        offset = 0; %case with no offset but change in Temperature
    case 3 
        params.a.T_eV = 1e3;
        offset = 10^7; %case with offset and no change in Temperature
end

maxwell = @(v,x,y) a.n/(pi^3/2.*y^3).*exp(-((v-x)/y).^2);

%% Define the dimensions
 
dim_v = DIMENSION(0,3e7);
dim_z = DIMENSION(0,pi);

dim_v.moment_dV = @(v,p,t) 2.*pi.*v.^2;
dim_z.moment_dV = @(z,p,t) sin(z);

dimensions = {dim_v, dim_z};
num_dims = numel(dimensions);

%% Define the analytic solution (optional)

soln1 = new_md_func(num_dims,{ ...    
    params.soln_v, ...
    params.init_cond_z, ...
    });
solutions = {soln1};

%% Define the initial conditions

ic1 = new_md_func(num_dims,{params.init_cond_v,params.init_cond_z,params.init_cond_t});
initial_conditions = {ic1};

%% Define the boundary conditions

BCL = new_md_func(num_dims,{...
    @(v,p,t) v.*0, ... 
    params.init_cond_z, ...
    params.init_cond_s});

BCR = ic1;

%% Define the terms of the PDE

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
termC_z = SD_TERM({pterm1,pterm2,pterm3,pterm4});
termC   = MD_TERM(num_dims,{termC_z,[]});

% term V1 == 1/v^2 d/dv(v^3(m_a/(m_a + m_b))nu_s f))
% term V1 == g(v) q(v)      [mass, g(v) = 1/v^2,  BC N/A]
% q(v) == d/dv(g2(v)f(v))   [grad, g2(v) = v^3(m_a/(m_a + m_b))nu_s, BCL= N, BCR=N]

g1 = @(v,p,t,dat) 1./v.^2;
g2 = @(v,p,t,dat) v.^3*p.a.m.*p.nu_s(v,p.a,p.b)./(p.a.m + p.b.m);

pterm1  = MASS(g1);
pterm2  = GRAD(num_dims,g2,-1,'N','D', BCL, BCR);
termV_s = SD_TERM({pterm1,pterm2});
termV1  = MD_TERM(num_dims,{termV_s,[]});

% term V2 == 1/v^2 d/dv(v^4*0.5*nu_par*d/dv(f))
% term V2 == g(v) q(v)      [mass, g(v) = 1/v^2,  BC N/A]
% q(v) == d/dv(g2(v)r(v))   [grad, g2(v) = v^4*0.5*nu_par, BCL= D, BCR=D]
% r(v) = d/dv(g3(v)f)       [grad, g3(v) = 1, BCL=N, BCR=N]

g1 = @(v,p,t,dat) 1./v.^2;
g2 = @(v,p,t,dat) v.^4.*0.5.*p.nu_par(v,p.a,p.b);
g3 = @(v,p,t,dat) v.*0 + 1;

pterm1      = MASS(g1);
pterm2      = GRAD(num_dims,g2,-1,'D','N');
pterm3      = GRAD(num_dims,g3,+1,'N','D', BCL, BCR);
termV_par   = SD_TERM({pterm1,pterm2,pterm3});
termV2      = MD_TERM(num_dims,{termV_par,[]});

terms = {termV1,termV2,termC};


%% Define sources 

sources = {};

%% Define function to set time step
    function dt=set_dt(pde,CFL)    
        Lmax = pde.dimensions{1}.max;
        LevX = pde.dimensions{1}.lev;
        dt = Lmax/2^LevX*CFL;
    end


%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,[],initial_conditions,solutions);

end
