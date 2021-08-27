function pde = mirror2_collision_div(opts)
% Two-dimensional magnetic mirror from the FP paper - evolution of t*he ion velocity dependence
% of f in the presence of Coulomb collisions with background electrons
% Also applied to the 2D test for CQL4D equations
%
% df/dt == sum_{b} ( nu_D/(2*sin(z)) d/dth ( sin(th) 1/v d (v f) /dth ) + 
%  1/v^2 (d/dv(v^2[(m_a/(m_a + m_b))v nu_s f) + 0.5*nu_par*v^2*d/dv(f)])
%
% v,th are spherical coordinates (v,th,phi), so the div and grad look like 
%
% div[] = 1/v^2 * d/dv * v^2[] + 1/sin(th) d/dth (sin(th)[]), grad[] =
% d/dv[],1/v * d/dth[]
%
% and the volument_element dV = v^2 sin(th)
%
%
% df/dt == %  sum_{b} (( 1/v^2 (d/dv (v^2[0.5*nu_par*v^2*d/dv(f)]) %%term1 
%                 + 1/v^2 (d/dv(v^2[v (m_a/(m_a + m_b))nu_s f)]) %%term2
%                 + nu_D/(2 sin(z)) d/dth ( sin(th) df/dth )) %%term3 
%
% split into three div(flux) terms (term1, term2, and term3)
%
% term1 is done using SLDG defining A(v)= sqrt(0.5*nu_par*v^2)
%
% eq1 :  df/dt == div(A(v) * q)        [pterm1: div (g(v)=A(v),+1, BCL=?, BCR=?)]
% eq2 :      q == A(v) * grad(f)       [pterm2: grad(g(v)=A(v),-1, BCL=D, BCR=N)]
%
%
% term2 is a div using B(v) = v (m_a/(m_a + m_b))nu_s
%
% eq1 :  df/dt == div(B(v) * f)       [pterm1: div(g(p)=B(v),+1, BCL=?, BCR=?]
%
%
% term3 is done using SLDG defining C(v) = sqrt(v*nu_D(v)/2)
%
% eq1 :  df/dt == div(C(v) * q)   [pterm1: div (g(v)=C(v),+1, BCL=D, BCR=D]
% eq2 :      q == C(v) * grad(f)  [pterm2: grad(g(v)=C(v),-1, BCL=N, BCR=N]
% Run with
%
% asgard(@mirror2_collision_div,'timestep_method','BE','case',3, 'lev', 4, 'deg', 3,'num_steps', 60, 'dt',5e-5, 'normalize_by_mass', true, 'calculate_mass', true)

params = mirror_parameters();

switch opts.case_
    case 1 
        params.a.T_eV = 0.05*params.b.T_eV; %Target temperature in Kelvin\
        offset = 10^6; %case with offset and change in Temperature
    case 2 
        params.a.T_eV = 0.05*params.b.T_eV;
        offset = 0; %case with no offset but change in Temperature
    case 3 
        params.a.T_eV = 10;
        params.a.E_eV = 3e3; %case with offset and no change in Temperature
end
maxwell = @(v,x,y) a.n/(pi^3/2.*y^3).*exp(-((v-x)/y).^2);

%% Define the dimensions
 
dim_v = DIMENSION(0,1e6);
dV_v = @(x,p,t,d) x.^2;
dim_v.moment_dV = dV_v;

dim_th = DIMENSION(0,pi);
dV_th = @(x,p,t,d) sin(x);
dim_th.moment_dV = dV_th;

dimensions = {dim_v, dim_th};
num_dims = numel(dimensions);

%% Define the analytic solution (optional)

soln1 = new_md_func(num_dims,{ ...    
    params.soln_v, ...
    params.soln_z, ...
    });
solutions = {soln1};

%% Define the initial conditions

ic1 = new_md_func(num_dims,{params.init_cond_v,params.init_cond_z,params.init_cond_t});
initial_conditions = {ic1};

%% Define the boundary conditions

BCL = new_md_func(num_dims,{...
    @(v,p,t) v.*0, ...
    params.boundary_cond_z, ...
    params.boundary_cond_t});

BCR = new_md_func(num_dims,{...
    params.boundary_cond_v, ...
    params.boundary_cond_z, ...
    params.boundary_cond_t});

%% Define the terms of the PDE

%%
% term1 is done using SLDG defining A(v)= sqrt(0.5*nu_par*v^2)
%
% eq1 :  df/dt == div(A(v) * q)        [pterm1: div (g(v)=A(v),+1, BCL=?, BCR=?)]
% eq2 :      q == A(v) * grad(f)       [pterm2: grad(g(v)=A(v),-1, BCL=D, BCR=N)]

dV_v = @(x,p,t,d) x.^2; 
dV_th = @(x,p,t,d) sin(x);

g1 = @(v,p,t,dat) v.*0 + 1;
pterm1 = MASS(g1,'','',dV_th);
term1_th = SD_TERM({pterm1,pterm1});

A = @(v,p) sqrt(0.5.*v.^2.*(p.nu_par(v,p.a,p.b) + p.nu_par(v,p.a,p.b2)));
g1 = @(v,p,t,dat) A(v,p);
g2 = @(v,p,t,dat) A(v,p);

pterm1 = DIV (num_dims,g1,'',+1,'D','N','','','',dV_v);
pterm2 = GRAD(num_dims,g2,'',-1,'N','D','','','',dV_v);
term1_v = SD_TERM({pterm1,pterm2});
term1   = MD_TERM(num_dims,{term1_v,term1_th});

% term2 is a div using B(v) = v (m_a/(m_a + m_b))nu_s
%
% eq1 :  df/dt == div(B(v) * f)       [pterm1: div(g(p)=B(v),+1, BCL=?, BCR=?]

B = @(v,p) v*p.a.m.*(p.nu_s(v,p.a,p.b)./(p.a.m + p.b.m) + p.nu_s(v,p.a,p.b2)./(p.a.m + p.b2.m));
g3 = @(v,p,t,dat) B(v,p);

pterm1 = DIV(num_dims,g3,'',-1,'N','D','','','',dV_v);

% term3 is done using SLDG defining C(v) = sqrt(nu_D(v)/2)
%
% eq1 :  df/dt == div( q)   [pterm1: div (g(v)=C(v),+1, BCL=D, BCR=D]
% eq2 :      q ==  grad(f)  [pterm2: grad(g(v)=C(v),-1, BCL=N, BCR=N]

term2_v = SD_TERM({pterm1});
term2   = MD_TERM(num_dims,{term2_v,[]});

dV_v = @(x,p,t,d) x.^2; %changing for MASS term
dV_th = @(x,p,t,d) sin(x);

C = @(v,p) sqrt((p.nu_D(v,p.a,p.b) + p.nu_D(v,p.a,p.b2))/2);
g4 = @(v,p,t,dat) C(v,p);
pterm1 = MASS(g4,[],[],dV_v);
term3_v = SD_TERM({pterm1,pterm1});

D = @(v,p) v.*0 + 1;
g5 = @(v,p,t,dat) D(v,p);
pterm1 = DIV (num_dims,g5,'',+1,'D','N','','','',dV_th);
pterm2 = GRAD(num_dims,g5,'',-1,'N','D','','','',dV_th);
term3_th = SD_TERM({pterm1,pterm2});
term3   = MD_TERM(num_dims,{term3_v,term3_th});

terms = {term1,term2,term3};


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