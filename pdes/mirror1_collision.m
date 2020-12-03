function pde = mirror1_collision(opts)

%Testing the collision operator for the 1D form of the mirror equations
%
% df/dt = 1/v^2 (flux_v)
%
% flux_v == v^3[(m_a/(m_a + m_b))nu_s f) + 0.5*nu_par*v*d/dv(f)]
%
% df/dt = 1/v^2 (d/dv(v^3[(m_a/(m_a + m_b))nu_s f) + 0.5*nu_par*v*d/dv(f)])*v^2
%
% Run with
%
% asgard(@mirror1_collision,'timestep_method','BE','dt',1e-5,'num_steps',50,'lev',5,'deg',2,'normalize_by_mass',false, 'calculate_mass', false)

params = mirror_parameters();

%% Define the dimensions

dim_v = DIMENSION(0,1e6);
dim_v.jacobian = @(v,p,t) 4*pi.*v.^2;
dimensions = {dim_v};
num_dims = numel(dimensions);

%% Define the analytic solution (optional)

%     function ret = solution(v,p,t)
%         ret =  p.soln_v(v,p,t);
%         if isfield(p,'norm_fac')
%             ret = p.norm_fac .* ret;
%         end
%     end
soln_v = @solution_v;
    function ret = solution_v(v,p,t)
        ret = p.maxwell(v,0,p.v_th(p.b.T_eV,p.a.m));
        if isfield(p,'norm_fac')
            ret = p.norm_fac .* ret;
        end
    end

soln1 = new_md_func(num_dims,{soln_v});
solutions = {soln1};

%% Define the initial conditions

ic_v = @(v,p,t) p.init_cond_v(v);
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

% term V1 == 1/v^2 d/dv(v^3(m_a/(m_a + m_b))nu_s f))
% term V1 == g(v) q(v)      [mass, g(v) = 1/v^2,  BC N/A]
% q(v) == d/dv(g2(v)f(v))   [grad, g2(v) = v^3(m_a/(m_a + m_b))nu_s, BCL= N, BCR=N]

g1 = @(v,p,t,dat) 1./v.^2;%v.*0 + 1; 
g2 = @(v,p,t,dat) v.^3*p.a.m.*p.nu_s(v,p.a,p.b)./(p.a.m + p.b.m);

pterm1 = MASS(g1);
pterm2  = GRAD(num_dims,g2,0,'N','D', BCL, BCR);
termV_s = SD_TERM({pterm1,pterm2});
termV1   = MD_TERM(num_dims,{termV_s});

%%
% term V2 == 1/v^2 d/dv(v^4*0.5*nu_par*d/dv(f))
% term V2 == g(v) q(v)      [mass, g(v) = 1/v^2,  BC N/A]
% q(v) == d/dv(g2(v)r(v))   [grad, g2(v) = v^4*0.5*nu_par, BCL= D, BCR=D]
% r(v) = d/dv(g3(v)f)       [grad, g3(v) = 1, BCL=N, BCR=N]

g1 = @(v,p,t,dat) 1./v.^2; %v.*0 + 1;
g2 = @(v,p,t,dat) v.^4.*0.5.*p.nu_par(v,p.a,p.b);
g3 = @(v,p,t,dat) v.*0 + 1;

pterm1 = MASS(g1);
pterm2 = GRAD(num_dims,g2,-1,'D','N');
pterm3 = GRAD(num_dims,g3,+1,'N','D', BCL, BCR);
termV_par = SD_TERM({pterm1,pterm2,pterm3});
termV2   = MD_TERM(num_dims,{termV_par});

terms = {termV1,termV2};

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