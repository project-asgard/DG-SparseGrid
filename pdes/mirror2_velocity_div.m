function pde = mirror2_velocity_div(opts)
% Two-dimensional magnetic mirror from the FP paper - evolution of the ion velocity dependence
% of f in the presence of Coulomb collisions with background electrons
% Also applied to the 2D test for CQL4D equations

% df/dt == -1/v^2 d(v^2 ZeE/m cos(th) f)/dv + 1/(vsin(th)) d(sin(th) ZeE/m sin(th) f)/dth 
%            + sum_{b} ( nu_D/(2*sin(th)) d/dth ( sin(th) 1/v d (v f) /dth ) + 
%             1/v^2 (d/dv(v^2[(m_a/(m_a + m_b))v nu_s f) + 0.5*nu_par*v^2*d/dv(f)])
%
% v,th are spherical coordinates (v,th,phi), so the div and grad look like 
%
% div[] = 1/v^2 * d/dv * v^2[] + 1/sin(th) d/dth (sin(th)[]), grad[] =
% d/dv[],1/v * d/dth[]
%
% and the volument_element dV = v^2 sin(th)
%

% df/dt ==      -1/v^2 d(v^2 ZeE/m cos(th) f)/dv %term1
%                + 1/(vsin(th)) d(sin(th) ZeE/m sin(th) f)/dth  %term2
%                + % sum_b (( 1/v^2 (d/dv (v^2[0.5*nu_par*v^2*d/dv(f)]) %%term3 
%                 + 1/v^2 (d/dv(v^2[v (m_a/(m_a + m_b))nu_s f)]) %%term4
%                 + nu_D/(2 sin(z)) d/dth ( sin(th) df/dth )) %%term5 
%
% split into five div terms (term1,term2, term3, term4, and term5)
%
% term1 is done combining mass and div defining F(th) = -ZeE/m cos(th) and
% G(v) = 1
%
% eq1 : df/dt == div(F(th) f)      [pterm1: div(g(v)=G(v),-1, BCL=D,BCR=N)]
%
% term2 is a simple div term, defining K(th) = ZeE/m sin(th)
%
% eq1 : df/dt == div(K(th) f)     [pterm1:div(g(th)=K(th),-1, BCL=N,BCR=N)]
%
% term3 is done combining mass and div defining C(v) = (nu_D(v)/2) 
%
% eq1 :  df/dt == div( q)   [pterm1: div (g(v)=D(v),+1, BCL=D, BCR=D]
% eq2 :      q ==  C(v) grad(f)  [pterm2: grad(g(v)=D(v),-1, BCL=N, BCR=N]
%
% term4 is a div using B(v) = v (m_a/(m_a + m_b))nu_s
%
% eq1 :  df/dt == div(B(v) * f)       [pterm1: div(g(p)=B(v),+1, BCL=?, BCR=?]
%
%
% term5 is done using SLDG defining C(v) = sqrt(v*nu_D(v)/2)
%
% eq1 :  df/dt == div(C(v) * q)   [pterm1: div (g(v)=C(v),+1, BCL=D, BCR=D]
% eq2 :      q == C(v) * grad(f)  [pterm2: grad(g(v)=C(v),-1, BCL=N, BCR=N]
% Run with
%
% asgard(@mirror2_velocity_div,'timestep_method','BE','case',3, 'grid_type','SG', 'lev', 3, 'deg', 3,'num_steps', 10, 'dt',5e-5, 'normalize_by_mass', true, 'calculate_mass', true,'save_output',true,'update_params_each_timestep', true)

% Run with
%
% asgard(@mirror2_velocity_div,'timestep_method','BE','case',3,'dt',1e-6)

params = mirror_parameters(opts);


%% Define the dimensions
 
dim_v = DIMENSION(0,10*params.a.vth);
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
    params.boundary_cond_v, ...
    params.boundary_cond_z, ... %params.boundary_cond_z, ...
    params.boundary_cond_t});

BCR = new_md_func(num_dims,{...
    params.boundary_cond_v, ...
    params.boundary_cond_z, ...
    params.boundary_cond_t});

%% Define the terms of the PDE

% termE1a is done combining mass and div defining F(th) = (ZeE/m cos(th)) and
% G(v) = 1
%
% eq1 : df/dt == div(F(th) f)      [pterm1: div(g(v)=G(v),-1, BCL=D,BCR=N)]

 dV_v = @(x,p,t,d) x.^2; 
 dV_th = @(x,p,t,d) sin(x);

F = @(x,p) -cos(x).*p.a.Z.*p.e.*p.E./p.a.m;
g1 = @(x,p,t,dat) F(x,p).*(x>pi/2);
pterm1 = MASS(g1,'','',dV_th);
termE1_th = SD_TERM({pterm1});

G = @(v,p) v.*0 + 1;
g2 = @(v,p,t,dat) G(v,p);
pterm1 = DIV(num_dims,g2,'',+1,'D','N',BCL,'','',dV_v);
termE1_v = SD_TERM({pterm1});
termE1a = MD_TERM(num_dims,{termE1_v,termE1_th});

%termE1b is the same form as term1 but accounting for the flow in the
%opposite direction

F = @(x,p) -cos(x).*p.a.Z.*p.e.*p.E./p.a.m;
g1 = @(x,p,t,dat) F(x,p).*(x<pi/2);
pterm1 = MASS(g1,'','',dV_th);
termE1_th = SD_TERM({pterm1});

G = @(v,p) v.*0 + 1;
g2 = @(v,p,t,dat) G(v,p);
pterm1 = DIV(num_dims,g2,'',-1,'N','D','',BCR,'',dV_v);
termE1_v = SD_TERM({pterm1});
termE1b = MD_TERM(num_dims,{termE1_v,termE1_th});

% termE2 is a simple div term, defining K(th) = ZeE/m sin(th)
%
% eq1 : df/dt == div(K(th) f)     [pterm1:div(g(th)=K(th),-1, BCL=N,BCR=N)]

% dV_v = @(x,p,t,d) x.^2; %changing for MASS term
% dV_th = @(x,p,t,d) sin(x);

dV_v = @(x,p,t,d) x; %changing for MASS term
dV_th = @(x,p,t,d) sin(x);

g1 = @(x,p,t,dat) x.*0+1;
pterm1   =  MASS(g1,'','',dV_v);
termE2_v = SD_TERM({pterm1});

K = @(x,p) p.a.Z.*p.e.*p.E.*sin(x)./p.a.m;
g2 = @(x,p,t,dat) K(x,p);

pterm1 = DIV(num_dims,g2,'',-1,'N','N','','','',dV_th);
termE2_th = SD_TERM({pterm1});
termE2 = MD_TERM(num_dims,{termE2_v,termE2_th});

% dV_v = @(x,p,t,d) x.^2;
% dV_th = @(x,p,t,d) sin(x);

%%
% termC1 is done using SLDG defining A(v)= sqrt(0.5*nu_par*v^2)
%
% eq1 :  df/dt == div(A(v) * q)        [pterm1: div (g(v)=A(v),+1, BCL=?, BCR=?)]
% eq2 :      q == A(v) * grad(f)       [pterm2: grad(g(v)=A(v),-1, BCL=D, BCR=N)]

dV_v = @(x,p,t,d) x.^2; 
dV_th = @(x,p,t,d) sin(x);

g1 = @(v,p,t,dat) v.*0 + 1;
pterm1 = MASS(g1,'','',dV_th);
termC1_th = SD_TERM({pterm1,pterm1});

A = @(v,p) v.*sqrt(0.5*(p.nu_par(v,p.a,p.b) + p.nu_par(v,p.a,p.b2)));
g1 = @(v,p,t,dat) A(v,p);
g2 = @(v,p,t,dat) A(v,p);

pterm1 = DIV (num_dims,g1,'',+1,'D','N','','','',dV_v);
pterm2 = GRAD(num_dims,g2,'',-1,'N','D','','','',dV_v);
termC1_v = SD_TERM({pterm1,pterm2});
termC1   = MD_TERM(num_dims,{termC1_v,termC1_th});

% termC2 is a div using B(v) = v (m_a/(m_a + m_b))nu_s
%
% eq1 :  df/dt == div(B(v) * f)       [pterm1: div(g(p)=B(v),+1, BCL=?, BCR=?]

dV_v = @(x,p,t,d) x.^2; 
dV_th = @(x,p,t,d) sin(x);

A = @(v,p) v.*0 + 1;
g1 = @(v,p,t,dat) A(v,p);

pterm1 = MASS(g1,'','',dV_th);

termC2_th = SD_TERM({pterm1});

B = @(v,p) v.*p.a.m.*(p.nu_s(v,p.a,p.b)./(p.a.m + p.b.m) + p.nu_s(v,p.a,p.b2)./(p.a.m + p.b2.m));
g2 = @(v,p,t,dat) B(v,p);

pterm2 = DIV(num_dims,g2,'',-1,'N','D','','','',dV_v);

term2_v = SD_TERM({pterm2});
termC2   = MD_TERM(num_dims,{term2_v,termC2_th});

% termC3 is done using SLDG defining C(v) = sqrt(nu_D(v)/2)
%
% eq1 :  df/dt == div( q)   [pterm1: div (g(v)=C(v),+1, BCL=D, BCR=D]
% eq2 :      q ==  grad(f)  [pterm2: grad(g(v)=C(v),-1, BCL=N, BCR=N]


dV_v = @(x,p,t,d) x; %changing for MASS term
dV_th = @(x,p,t,d) sin(x);

C = @(v,p) v.*sqrt((p.nu_D(v,p.a,p.b) + p.nu_D(v,p.a,p.b2))/2);
g4 = @(v,p,t,dat) C(v,p);
pterm1 = MASS(g4,'','',dV_v);
termC3_v = SD_TERM({pterm1,pterm1});

D = @(v,p) v.*0 + 1;
g5 = @(v,p,t,dat) D(v,p);
pterm1 = DIV (num_dims,g5,'',+1,'D','D','','','',dV_th);
pterm2 = GRAD(num_dims,g5,'',-1,'N','N','','','',dV_th);
termC3_th = SD_TERM({pterm1,pterm2});
termC3   = MD_TERM(num_dims,{termC3_v,termC3_th});

terms = {termE1a,termE1b,termE2,termC1,termC2,termC3};

%% Define sources 

    function res = my_alpha(x,p,t)
%         disp(num2str(p.alpha_z(x)));
        res = p.alpha_z(x);
    end

sources = {};

switch opts.case_
    case 3
        source1_v = @(x,p,t,d) x.*0 + 1;%p.f0_v(x);
        source1_z = @(x,p,t,d) my_alpha(x,p,t);
        source1 = new_md_func(num_dims,{source1_v,source1_z});
        sources = {source1};
end
%% Define function to set time step
    function dt=set_dt(pde,CFL)    
        Lmax = pde.dimensions{1}.max;
        LevX = pde.dimensions{1}.lev;
        dt = Lmax/2^LevX*CFL;
    end


%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,[],initial_conditions,solutions);

end
