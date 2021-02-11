function pde = mirror1_velocity(opts)
% One-dimensional magnetic mirror from the FP paper - evolution of the ion velocity dependence
% of f in the presence of Coulomb collisions with background electrons
% 
% df/dt v^2 == (v^2/2)df/dv v^2  - v*f*v^2 + 1/v^2 (d/dv(flux_v))*v^2
%
% flux_v == v^3[(m_a/(m_a + m_b))nu_s f) + 0.5*nu_par*v*d/dv(f)]
%
% df/dt == (v^2)/2 )df/dv*v^2 -v*f*v^2 + 1/v^2 (d/dv(v^3[(m_a/(m_a +
% m_b))nu_s f) + 0.5*nu_par*v*d/dv(f)])*v^2
%
%
% Run with
%
% asgard(@mirror1_velocity,'timestep_method','BE','dt',1e-7,'num_steps',50,'lev',5,'deg',2,'normalize_by_mass',false, 'calculate_mass', false)

params = mirror_parameters();

%% Define the dimensions

dim_v = DIMENSION(0,1e6);
dim_v.jacobian = @(v,p,t) v.^2;
dimensions = {dim_v};
num_dims = numel(dimensions);

%% Define the analytic solution (optional)

    function ret = solution(v,p,t)
        ret =  p.soln_v(v,p,t);
        if isfield(p,'norm_fac')
            ret = p.norm_fac .* ret;
        end
    end

soln1 = new_md_func(num_dims,{@solution});
solutions = {soln1};

%% Define the initial conditions

ic_v = @(v,p,t) p.a.n.*p.init_cond_v(v);
ic1 = new_md_func(num_dims,{ic_v});
initial_conditions = {ic1};

%% Define the boundary conditions

%% Define the boundary conditions

BCL = new_md_func(num_dims,{...
    @(v,p,t) v.*0, ...
    params.boundary_cond_t});
 

BCR = new_md_func(num_dims,{...
    params.boundary_cond_v, ... 
    params.boundary_cond_t});

%% Define the terms of the PDE

%% 
% -E*Z_a/m_a d/dz(f)

%g1 = @(v,p,t,dat) -E.*Z_a/m_a;
%pterm1  = GRAD(num_dims,g1,-1,'N','N');
%termE_v = SD_TERM({pterm1});
%termE   = MD_TERM(num_dims,{termE_v});

%% 
% LHS_term == df/dt v^2
g1 = @(v,p,t,dat) dim_v.jacobian(v,p,t);
pterm1 = MASS(g1);
LHS_term_z = SD_TERM({pterm1});
LHS_term = MD_TERM(num_dims,{LHS_term_z});
LHS_terms = {LHS_term};


% term V1 == v^2/2 df/dv*v^2
% term V1 == g(v) q(v) 
% g(v) = g1(v) [mass, g1(z) = v^2/2, BC = N/A]
% q(v) = d/dv (g2(v) f) [grad, g2(v) = 1, BCL=N, BCR=D ]
g1 = @(v,p,t,dat) v.^2.*dim_v.jacobian(v,p,t)./2;
g2 = @(v,p,t,dat) v.*0 + 1;
pterm1 = MASS(g1);
pterm2 = GRAD(num_dims,g2,-1,'N','D', BCL, BCR);
termV_v = SD_TERM({pterm1,pterm2});
termV1 = MD_TERM(num_dims,{termV_v});

%term V2 == -v*f*v^2
%term V2 == g1(v)*f [mass, g1(v) = -v, BC = N/A]
g1 = @(v,p,t,dat) -v.*dim_v.jacobian(v,p,t);
pterm1 = MASS(g1);
termV_v = SD_TERM({pterm1});
termV2 = MD_TERM(num_dims,{termV_v});

% term V3 == 1/v^2 d/dv(v^3(m_a/(m_a + m_b))nu_s f))*v.^2
% term V3 == g(v) q(v)      [mass, g(v) = 1/v^2,  BC N/A]
% q(v) == d/dv(g2(v)f(v))   [grad, g2(v) = v^3(m_a/(m_a + m_b))nu_s, BCL= N, BCR=N]

g1 = @(v,p,t,dat) dim_v.jacobian(v,p,t)./v.^2;
g2 = @(v,p,t,dat) v.^3*p.a.m.*p.nu_s(v,p.a,p.b)./(p.a.m + p.b.m);

pterm1 = MASS(g1);
pterm2  = GRAD(num_dims,g2,0,'N','D', BCL, BCR);
termV_s = SD_TERM({pterm1,pterm2});
termV3   = MD_TERM(num_dims,{termV_s});

%%
% term V4 == 1/v^2 d/dv(v^4*0.5*nu_par*d/dv(f))*v^2
% term V4 == g(v) q(v)      [mass, g(v) = 1/v^2,  BC N/A]
% q(v) == d/dv(g2(v)r(v))   [grad, g2(v) = v^4*0.5*nu_par, BCL= D, BCR=D]
% r(v) = d/dv(g3(v)f)       [grad, g3(v) = 1, BCL=N, BCR=N]

g1 = @(v,p,t,dat) dim_v.jacobian(v,p,t)./v.^2;
g2 = @(v,p,t,dat) v.^4.*0.5.*p.nu_par(v,p.a,p.b);
g3 = @(v,p,t,dat) v.*0 + 1;

pterm1 = MASS(g1);
pterm2 = GRAD(num_dims,g2,+1,'D','N');
pterm3 = GRAD(num_dims,g3,-1,'N','D', BCL, BCR);
termV_par = SD_TERM({pterm1,pterm2,pterm3});
termV4   = MD_TERM(num_dims,{termV_par});

terms = {termV1,termV2,termV3,termV4};

%% Define sources 

sources = {};

%% Define function to set time step
    function dt=set_dt(pde,CFL)       
        Lmax = pde.dimensions{1}.max;
        LevX = pde.dimensions{1}.lev;
        dt = Lmax/2^LevX*CFL;
    end

%% Construct PDE

pde = PDE(opts,dimensions,terms,LHS_terms,sources,params,@set_dt,[],initial_conditions,solutions);

end
