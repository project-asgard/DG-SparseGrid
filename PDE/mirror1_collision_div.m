function pde = mirror1_collision_div(opts)

%Testing the collision operator for the 1D form of the mirror equations
%
% df/dt == 1/v^2 (d/dv(v^2[v (m_a/(m_a + m_b))nu_s f) + 0.5*nu_par*v^2*d/dv(f)])
%       == 1/v^2 (d/dv (v^2[0.5*nu_par*v^2*d/dv(f)]) %%term1
%           + 1/v^2 (d/dv(v^2[v (m_a/(m_a + m_b))nu_s f)]) %%term2
       %
% v is a spherical coordinate (v,th,phi), so the div and grad look like 
%
% div[] = 1/v^2 * d/dv * v^2[], grad[] = d/dv[]
%
% and the volument_element dV = v^2
%
% 
% split into two div(flux) terms (term1 and term2)
%
% term1 is done using SLDG defining A(v)= sqrt(0.5*nu_par*v^2)
%
% eq1 :  df/dt == div(A(v) * q)        [pterm1: div (g(v)=A(v),+1, BCL=?, BCR=?)]
% eq2 :      q == A(v) * grad(f)       [pterm2: grad(g(v)=A(v),-1, BCL=D, BCR=N)]
%
% term2 is a div using B(v) = v (m_a/(m_a + m_b))nu_s
%
% eq1 :  df/dt == div(B(v) * f)       [pterm1: div(g(p)=B(v),+1, BCL=?, BCR=?]
%
%
% Run with
%
% asgard(@mirror1_collision_div,'timestep_method','BE','dt',1e-5,'num_steps',50,'lev',5,'deg',2,'normalize_by_mass',false, 'calculate_mass', false)

params = mirror_parameters();

%% Define the dimensions

dim_v = DIMENSION(0,1e6);
dV_v = @(x,p,t,d) x.^2;
dim_v.moment_dV = dV_v;
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

%% LHS terms (mass only)

LHS_terms = {};

%% Define the boundary conditions

%% Define the boundary conditions

BCL = new_md_func(num_dims,{...
    @(v,p,t) v.*0, ...
    params.boundary_cond_t});
 

BCR = new_md_func(num_dims,{...
    params.boundary_cond_v, ... 
    params.boundary_cond_t});


%%
% term1 is done using SLDG defining A(v)= sqrt(0.5*nu_par*v^2)
%
% eq1 :  df/dt == div(A(v) * q)        [pterm1: div (g(v)=A(v),+1, BCL=?, BCR=?)]
% eq2 :      q == A(v) * grad(f)       [pterm2: grad(g(v)=A(v),-1, BCL=D, BCR=N)]

A = @(v,p) sqrt(0.5.*v.^2.*p.nu_par(v,p.a,p.b));
g1 = @(v,p,t,dat) A(v,p);
g2 = @(v,p,t,dat) A(v,p);

pterm1 = DIV (num_dims,g1,'',+1,'D','N','','','',dV_v);
pterm2 = GRAD(num_dims,g2,'',-1,'N','D','','','',dV_v);
term1_v = SD_TERM({pterm1,pterm2});
term1   = MD_TERM(num_dims,{term1_v});

% term2 is a div using B(v) = v (m_a/(m_a + m_b))nu_s
%
% eq1 :  df/dt == div(B(v) * f)       [pterm1: div(g(p)=B(v),+1, BCL=?, BCR=?]

B = @(v,p) v*p.a.m.*p.nu_s(v,p.a,p.b)./(p.a.m + p.b.m);
g3 = @(v,p,t,dat) B(v,p);

pterm3 = DIV(num_dims,g3,'',-1,'N','D','','','',dV_v);

term2_v = SD_TERM({pterm3});
term2   = MD_TERM(num_dims,{term2_v});

terms = {term1,term2};

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