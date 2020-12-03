function pde = mirror3(opts)
% Three-dimensional magnetic mirror from the FP paper - evolution of the ion velocity and spatial dependence
% of f in the presence of Coulomb collisions with background electrons
% 
% df/dt == -dB/ds f - vcos(z)df/dz + nu_D/(2*sin(z)) d/dz ( sin(z) df/dz ) + 1/v^2 (d/dv(flux_v))
%a
% flux_v == v^3[(m_a/(m_a + m_b))nu_s f) + 0.5*nu_par*v*d/dv(f)]
%
% Run with
%
% asgard(@mirror3,'timestep_method','BE','case',3)

% choose a case
%
% case 1 - maxwellian offset and different Temp.
% case 2 - max. no offset and different Temp.
% case 3 - max. offset and same Temp.

params = mirror_parameters();

switch opts.case_
    case 1
        params.a.T_eV = 0.05*params.b.T_eV; %Target temperature in Kelvin
        params.init_cond_v = @(v,p,t) params.maxwell(v,params.v_th(params.a.T_eV,params.a.m),1e6);
    case 2
        params.a.T_eV = 250*params.b.T_eV;
        params.init_cond_v = @(v,p,t) params.maxwell(v,0, 1e6);
    case 3
        params.B_func = params.B_func2; %setting magnetic field to loop function
        params.dB_ds = params.dB_ds2;
        params.init_cond_s = @(s,p,t) exp(s);
end

%% Define the dimensions

dim_v = DIMENSION(0,2e7);
dim_z = DIMENSION(0,pi/2);
dim_s = DIMENSION(0,2);

dim_v.jacobian = @(v,p,t) 2.*pi.*v.^2;
dim_z.jacobian = @(z,p,t) sin(z);
dim_s.jacobian = @(s,p,t) s.*0 + 1;

dimensions = {dim_v,dim_z,dim_s};
num_dims = numel(dimensions);

%% Define the analytic solution (optional)

soln1 = new_md_func(num_dims,{params.soln_v,params.soln_z,params.soln_s});
solutions = {soln1};

%% Define the initial conditions

ic1 = new_md_func(num_dims,{params.init_cond_v,params.init_cond_z,params.init_cond_s});
initial_conditions = {ic1};

%% Define the boundary conditions

BCL = new_md_func(num_dims,{...
    @(v,p,t) v.*0, ... 
    params.boundary_cond_z, ...
    params.boundary_cond_s, ...
    params.boundary_cond_t});

BCR = new_md_func(num_dims,{...
    params.boundary_cond_v, ... 
    params.boundary_cond_z, ...
    params.boundary_cond_s, ...
    params.boundary_cond_t});

%% Define the terms of the PDE

%% Advection  term
% termS1 == -vcos(z)*df/ds
% termS1 == q(v)*r(z)*w(s)
% q(v) == g1(v)  [mass, g1(p) = v,  BC N/A]
% r(z) == g2(z) [mass, g2(z) = cos(z),  BC N/A]
% w(s) == d/ds g3(s) f [grad, g3(s) = -1, BCL= D, BCR=D]
g1 = @(v,p,t,dat) -v;
g2 = @(z,p,t,dat) cos(z);
g3 = @(s,p,t,dat) s.*0 + 1;

pterm1 = MASS(g1);
pterm2 = MASS(g2);
pterm3 = GRAD(num_dims,g3,-1,'D','D', BCL, BCR);

term1_v = SD_TERM({pterm1});
term1_z = SD_TERM({pterm2});
term1_s = SD_TERM({pterm3});
termS1  = MD_TERM(num_dims,{term1_v,term1_z,term1_s});

%% Mass term
% termS2 == -vcos(z)dB/ds f
% termS1 == q(v)*r(z)*w(s)
% q(v) == g1(v)  [mass, g1(p) = v,  BC N/A]
% r(z) == g2(z) [mass, g2(z) = cos(z),  BC N/A]
% w(s) == g3(s) f [mass, g3(s) = -dB/ds/B, BCL= D, BCR=D]

g1 = @(v,p,t,dat) -v;
g2 = @(z,p,t,dat) cos(z);
g3 = @(s,p,t,dat) p.dB_ds(s)./p.B_func(s);

pterm1 = MASS(g1);
pterm2 = MASS(g2);
pterm3 = MASS(g3);

termB1_v = SD_TERM({pterm1});
termB1_z = SD_TERM({pterm2});
termB1_s = SD_TERM({pterm3});

termS2 = MD_TERM(num_dims,{termB1_v,termB1_z,termB1_s});

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
pterm3  = GRAD(num_dims,g3,+1,'D','N');
pterm4  = GRAD(num_dims,g4,-1,'N', 'D', BCL, BCR);
termC_v = SD_TERM({pterm1});
termC_z = SD_TERM({pterm2,pterm3,pterm4});
termC   = MD_TERM(num_dims,{termC_v,termC_z,[]});

% term V1 == 1/v^2 d/dv(v^3(m_a/(m_a + m_b))nu_s f))
% term V1 == g(v) q(v)      [mass, g(v) = 1/v^2,  BC N/A]
% q(v) == d/dv(g2(v)f(v))   [grad, g2(v) = v^3(m_a/(m_a + m_b))nu_s, BCL= N, BCR=D]

g1 = @(v,p,t,dat) 1./v.^2;
g2 = @(v,p,t,dat) v.^3*p.a.m.*p.nu_s(v,p.a,p.b)./(p.a.m + p.b.m);

pterm1  = MASS(g1);
pterm2  = GRAD(num_dims,g2,-1,'N','D', BCL, BCR);
termV_v = SD_TERM({pterm1,pterm2});
termV1  = MD_TERM(num_dims,{termV_v,[],[]});

%%
% term V2 == 1/v^2 d/dv(v^4*0.5*nu_par*d/dv(f))
% term V2 == g(v) q(v)      [mass, g(v) = 1/v^2,  BC N/A]
% q(v) == d/dv(g2(v)r(v))   [grad, g2(v) = v^4*0.5*nu_par, BCL= D, BCR=D]
% r(v) = d/dv(g3(v)f)       [grad, g3(v) = 1, BCL=N, BCR=N]

g1 = @(v,p,t,dat) 1./v.^2;
g2 = @(v,p,t,dat) v.^4.*0.5.*p.nu_par(v,p.a,p.b);
g3 = @(v,p,t,dat) v.*0 + 1;

pterm1      = MASS(g1);
pterm2      = GRAD(num_dims,g2,+1,'D','N');
pterm3      = GRAD(num_dims,g3,-1,'N','D', BCL, BCR);
termV_v     = SD_TERM({pterm1,pterm2,pterm3});
termV2      = MD_TERM(num_dims,{termV_v,[],[]});
terms = {termV1, termV2, termC, termS1, termS2};

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
