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

params = mirror_parameters();

switch opts.case_
    case 1 
        params.a.T_eV = 50;
        params.a.vth = params.v_th(params.a.T_eV,params.a.m);
        params.a.E_eV = 3e3; %case with offset and no change in Temperature
        params.init_cond_v = @(v,p,t) params.maxwell(v,params.a.v_beam,params.a.vth);
    case 2 
        m_e_cgs = 9.109*10^-28; %electron mass in g
        m_D_cgs = 3.3443*10^-24; %Deuterium mass in g
        m_He_cgs = 6.7*10^-24; %helium 4 mass in g 
        m_B_cgs = 1.82*10^-23; %Boron 11 mass in g
        temp_cgs = 1.6022e-10; %temperature in erg
%         params_cgs.a.vth = sqrt(2*temp_cgs/m_e_cgs);
%         params_cgs.b.vth = sqrt(2*temp_cgs/m_e_cgs);
        params_cgs.a.m = m_e_cgs; %beam is electrons
        params_cgs.b.m = m_e_cgs; %background is electrons
        params_cgs.b2.m = m_D_cgs;
        params_cgs.a.Z = -1;
        params_cgs.b.Z = -1;
        params_cgs.b2.Z = 1;
        params_cgs.e = 4.803*10^-10; %charge in Fr
        params_cgs.E = 2.6e-5; %E field in statvolt/cm
%         params.a.vth = 0.01*params_cgs.a.vth; %converting to m/s
%         params.b.vth = 0.01*params_cgs.b.vth;
        params.a.m = 0.001*params_cgs.a.m; %converting to kg
        params.b.m = 0.001*params_cgs.b.m; 
        params.b2.m = 0.001*params_cgs.b2.m;
        params.a.Z = params_cgs.a.Z;
        params.b.Z = params_cgs.b.Z;
        params.b2.Z = params_cgs.b2.Z;
        params.a.n = 5e19;
        params.a.T_eV = 5111;
        params.ln_delt = 20;
        params.a.vth = params.v_th(params.a.T_eV,params.a.m)/sqrt(2);
        params.b.vth = params.v_th(params.b.T_eV,params.b.m)/sqrt(2);
        E_dreicer_si = params.a.n.*params.e^3*params.ln_delt/(4*pi*params.eps0^2*params.a.m ... 
            *params.a.vth^2);
        params.E = -10^-4*E_dreicer_si;
        %vel_norm = @(v,vth) v./vth; %normalized velocity to thermal velocity
        params.maxwell = @(v,offset,vth) params.a.n/(pi.^(3/2)*vth^3).*exp(-((v-offset)/vth).^2);
        params.init_cond_v = @(v,p,t) params.maxwell(v,0,params.a.vth);
        %params_cgs.nu_ab0  = @(a,b) b.n * params_cgs.e^4 * a.Z^2 * b.Z^2 * params_cgs.ln_delt / (pi^3/2.*a.m^2*b.vth^3); %scaling coefficient
    case 3 
        n_cgs = 8e14; %equilibrium density in cm.^-3
        m_e_cgs = 9.109*10^-28; %electron mass in g
        m_D_cgs = 3.3443*10^-24; %Deuterium mass in g
        m_He_cgs = 6.7*10^-24; %helium 4 mass in g 
        m_B_cgs = 1.82*10^-23; %Boron 11 mass in g
        m_Ne_cgs = 3.3509177*10^-23; %Neon mass in g
        temp_cgs = 1.6022e-10; %temperature in erg
        params_cgs.a.n = n_cgs;
        params_cgs.b.n = n_cgs;
        params_cgs.b2.n = n_cgs;
%         params_cgs.a.vth = sqrt(2*temp_cgs/m_e_cgs);
%         params_cgs.b.vth = sqrt(2*temp_cgs/m_e_cgs);
        params_cgs.a.m = m_e_cgs; %beam is electrons
        params_cgs.b.m = m_D_cgs; %background ions
        params_cgs.b2.m = m_e_cgs; %background electrons
        params_cgs.a.Z = -1;
        params_cgs.b.Z = 1;
        params_cgs.b2.Z = -1;
        params_cgs.e = 4.803*10^-10; %charge in Fr
        params_cgs.E = 2.6e-5; %E field in statvolt/cm
        params.a.n = 10^6*params_cgs.a.n;%converting to m^-3
        params.b.n = 10^6*params_cgs.b.n;
        params.b2.n = 10^6*params_cgs.b2.n;
%         params.a.vth = 0.01*params_cgs.a.vth; %converting to m/s
%         params.b.vth = 0.01*params_cgs.b.vth;
        params.a.m = 0.001*params_cgs.a.m; %converting to kg
        params.b.m = 0.001*params_cgs.b.m; 
        params.b2.m = 0.001*params_cgs.b2.m;
        params.a.Z = params_cgs.a.Z;
        params.b.Z = params_cgs.b.Z;
        params.b2.Z = params_cgs.b2.Z;
        %params.E = 2.9979*10^4*params_cgs.E; %converting to V/m
        params.a.E_eV = 1000;
        params.a.T_eV = 5.11*10^3;
        params.b.T_eV = params.a.T_eV;
        params.b2.T_eV = params.a.T_eV;
        params.a.vth = params.v_th(params.a.T_eV,params.a.m);
        params.b.vth = params.v_th(params.b.T_eV,params.b.m);
        params.b2.vth = params.v_th(params.b2.T_eV,params.b2.m);
        params.ln_delt = 15;
        E_dreicer_si = params.a.n.*params.e^3*params.ln_delt/(2*pi*params.eps0^2*params.a.m ... 
            *params.a.vth^2);
        frac = 1e-6;
        params.E = frac*E_dreicer_si;
        %vel_norm = @(v,vth) v./vth; %normalized velocity to thermal velocity
        params.maxwell = @(v,offset,vth) params.a.n/(pi.^(3/2)*vth^3).*exp(-((v-offset)/vth).^2);
        params.soln_v = @(v,p,t) solution_v(v,p,t);
        params.f0_v = @(v) params.maxwell(v,0,params.a.vth);
        params.init_cond_v = @(v,p,t) params.f0_v(v);
        %params_cgs.nu_ab0  = @(a,b) b.n * params_cgs.e^4 * a.Z^2 * b.Z^2 * params_cgs.ln_delt / (pi^3/2.*a.m^2*b.vth^3); %scaling coefficient
        %params.eps0 = 1/(4*pi);
end
    function ret = phi(x)
        ret = erf(x);
    end

    function ret = dphidx(x)
        ret = 2./sqrt(pi) * exp(-x.^2);
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

    function ret = solution_v(v,p,t)
        ret = params.a.n/(pi^3/2.*params.v_th(params.b.T_eV,params.a.m).^3).*...
            exp(-(v./params.v_th(params.b.T_eV,params.a.m)).^2);
        if isfield(p,'norm_fac')
            ret = p.norm_fac .* ret;
        end
    end
maxwell = @(v,x,y) a.n/(pi^3/2.*y^3).*exp(-((v-x)/y).^2);

%% Define the dimensions
 
dim_v = DIMENSION(0,6*params.a.vth);
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
    @(v,p,t) v.*0, ...
    params.boundary_cond_z, ...
    params.boundary_cond_t});

%% Define the terms of the PDE

% termE1a is done combining mass and div defining F(th) = (ZeE/m cos(th)) and
% G(v) = 1
%
% eq1 : df/dt == div(F(th) f)      [pterm1: div(g(v)=G(v),-1, BCL=D,BCR=N)]

 dV_v = @(x,p,t,d) x.^2; 
 dV_th = @(x,p,t,d) sin(x);

F = @(x,p) -cos(x).*params.a.Z.*params.e.*params.E./params.a.m;
g1 = @(x,p,t,dat) F(x,p).*(x>pi/2);
pterm1 = MASS(g1,'','',dV_th);
termE1_th = SD_TERM({pterm1});

G = @(v,p) v.*0 + 1;
g2 = @(v,p,t,dat) G(v,p);
pterm1 = DIV(num_dims,g2,'',-1,'N','D','',BCR,'',dV_v);
termE1_v = SD_TERM({pterm1});
termE1a = MD_TERM(num_dims,{termE1_v,termE1_th});

%termE1b is the same form as term1 but accounting for the flow in the
%opposite direction

F = @(x,p) -cos(x).*params.a.Z.*params.e.*params.E./params.a.m;
g1 = @(x,p,t,dat) F(x,p).*(x<pi/2);
pterm1 = MASS(g1,'','',dV_th);
termE1_th = SD_TERM({pterm1});

G = @(v,p) v.*0 + 1;
g2 = @(v,p,t,dat) G(v,p);
pterm1 = DIV(num_dims,g2,'',+1,'D','N','',BCR,'',dV_v);
termE1_v = SD_TERM({pterm1});
termE1b = MD_TERM(num_dims,{termE1_v,termE1_th});

% termE2 is a simple div term, defining K(th) = ZeE/m sin(th)
%
% eq1 : df/dt == div(K(th) f)     [pterm1:div(g(th)=K(th),-1, BCL=N,BCR=N)]

% dV_v = @(x,p,t,d) x.^2; %changing for MASS term
% dV_th = @(x,p,t,d) sin(x);

dV_v = @(x,p,t,d) x; %changing for MASS term
dV_th = @(x,p,t,d) sin(x);

g1 = @(x,p,t,dat) 0*x+1;
pterm1   =  MASS(g1,'','',dV_v);
termE2_v = SD_TERM({pterm1});

K = @(x,p) params.a.Z.*params.e.*params.E.*sin(x)./params.a.m;
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