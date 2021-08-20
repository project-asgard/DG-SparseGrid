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
% asgard(@mirror2_velocity_div,'timestep_method','matrix_exponential','case',3, 'lev', 4, 'deg', 3,'num_steps', 60, 'dt',5e-5, 'normalize_by_mass', true, 'calculate_mass', true)

% Run with
%
% asgard(@mirror2_velocity_div,'timestep_method','BE','case',3,'dt',1e-6)

params_si = mirror_parameters();

switch opts.case_
    case 1 
        params_si.a.T_eV = 0.05*params_si.b.T_eV; %Target temperature in Kelvin\
        offset = 10^6; %case with offset and change in Temperature
    case 2 
        params_si.a.T_eV = 0.05*params_si.b.T_eV;
        offset = 0; %case with no offset but change in Temperature
    case 3 
        n_cgs = 8e14; %equilibrium density in cm.^-3
        m_e_cgs = 9.109*10^-28; %electron mass in g
        m_D_cgs = 3.3443*10^-24; %Deuterium mass in g
        m_He_cgs = 6.7*10^-24; %helium 4 mass in g 
        m_B_cgs = 1.82*10^-23; %Boron 11 mass in g
        temp_cgs = 1.6022e-10; %temperature in erg
        params_cgs.a.n = n_cgs;
        params_cgs.b.n = n_cgs;
        params_cgs.b2.n = n_cgs;
%         params_cgs.a.vth = sqrt(2*temp_cgs/m_e_cgs);
%         params_cgs.b.vth = sqrt(2*temp_cgs/m_e_cgs);
        params_cgs.a.m = m_e_cgs; %beam is electrons
        params_cgs.b.m = m_e_cgs; %background is electrons
        params_cgs.b2.m = m_B_cgs;
        params_cgs.a.Z = -1;
        params_cgs.b.Z = -1;
        params_cgs.b2.Z = 5;
        params_cgs.e = 4.803*10^-10; %charge in Fr
        params_cgs.E = 2.6e-5; %E field in statvolt/cm
        params_si.a.n = 10^6*params_cgs.a.n;%converting to m^-3
        params_si.b.n = 10^6*params_cgs.b.n;
        params_si.b2.n = 10^6*params_cgs.b2.n;
%         params_si.a.vth = 0.01*params_cgs.a.vth; %converting to m/s
%         params_si.b.vth = 0.01*params_cgs.b.vth;
        params_si.a.m = 0.001*params_cgs.a.m; %converting to kg
        params_si.b.m = 0.001*params_cgs.b.m; 
        params_si.b2.m = 0.001*params_cgs.b2.m;
        params_si.a.Z = params_cgs.a.Z;
        params_si.b.Z = params_cgs.b.Z;
        params_si.b2.Z = params_cgs.b2.Z;
        %params_si.E = 2.9979*10^4*params_cgs.E; %converting to V/m
        params_si.a.E_eV = 7.665;
        params_si.a.T_eV = 2/3*params_si.a.E_eV;
        params_si.b.T_eV = params_si.a.T_eV;
        params_si.b2.T_eV = params_si.a.T_eV;
        params_si.a.vth = params_si.v_th(params_si.a.T_eV,params_si.a.m);
        params_si.b.vth = params_si.v_th(params_si.b.T_eV,params_si.b.m);
        params_si.b2.v.th = params_si.v_th(params_si.b2.T_eV,params_si.b2.m);
        params_si.ln_delt = 15;
        E_dreicer_si = params_si.a.n.*params_si.e^3*params_si.ln_delt/(4*pi*params_si.eps0^2*params_si.a.m ... 
            *params_si.a.vth^2);
        params_si.E = 10^-6*E_dreicer_si;
        %vel_norm = @(v,vth) v./vth; %normalized velocity to thermal velocity
        params_si.maxwell = @(v,offset,vth) params_si.a.n/(pi.^(3/2)*vth^3).*exp(-((v-offset)/vth).^2);
        params_si.init_cond_v = @(v,p,t) params_si.maxwell(v,0,params_si.a.vth);
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

maxwell = @(v,x,y) a.n/(pi^3/2.*y^3).*exp(-((v-x)/y).^2);

%% Define the dimensions
 
dim_v = DIMENSION(0,10*params_si.a.vth);
dV_v = @(x,p,t,d) x.^2;
dim_v.moment_dV = dV_v;

dim_th = DIMENSION(0,pi);
dV_th = @(x,p,t,d) sin(x);
dim_th.moment_dV = dV_th;

dimensions = {dim_v, dim_th};
num_dims = numel(dimensions);

%% Define the analytic solution (optional)

soln1 = new_md_func(num_dims,{ ...    
    params_si.soln_v, ...
    params_si.soln_z, ...
    });
solutions = {soln1};

%% Define the initial conditions

ic1 = new_md_func(num_dims,{params_si.init_cond_v,params_si.init_cond_z,params_si.init_cond_t});
initial_conditions = {ic1};

%% Define the boundary conditions

BCL = new_md_func(num_dims,{...
    params_si.init_cond_v, ...
    params_si.boundary_cond_z, ...
    params_si.boundary_cond_t});

BCR = new_md_func(num_dims,{...
    @(v,p,t) v.*0, ...
    params_si.boundary_cond_z, ...
    params_si.boundary_cond_t});

%% Define the terms of the PDE

% term1 is done combining mass and div defining F(th) = (ZeE/m cos(th)) and
% G(v) = 1
%
% eq1 : df/dt == div(F(th) f)      [pterm1: div(g(v)=G(v),-1, BCL=D,BCR=N)]

% dV_v = @(x,p,t,d) x; %changing for MASS term
% dV_th = @(x,p,t,d) sin(x);

F = @(x,p) cos(x);
g1 = @(x,p,t,dat) F(x,p).*(x>pi/2);
pterm1 = MASS(g1,[],[],dV_th);
term1_th = SD_TERM({pterm1});

G = @(v,p) v.*0 - params_si.a.Z.*params_si.e.*params_si.E./params_si.a.m;
g2 = @(v,p,t,dat) G(v,p);
pterm1 = DIV(num_dims,g2,'',+1,'D','N',BCL,'','',dV_v);
term1_v = SD_TERM({pterm1});
term1a = MD_TERM(num_dims,{term1_v,term1_th});

%term1b is the same form as term1 but accounting for the flow in the
%opposite direction

F = @(x,p) cos(x);
g1 = @(x,p,t,dat) F(x,p).*(x<pi/2);
pterm1 = MASS(g1,[],[],dV_th);
term1_th = SD_TERM({pterm1});

G = @(v,p) v.*0 + params_si.a.Z.*params_si.e.*params_si.E./params_si.a.m;
g2 = @(v,p,t,dat) G(v,p);
pterm1 = DIV(num_dims,g2,'',-1,'N','D','',BCR,'',dV_v);
term1_v = SD_TERM({pterm1});
term1b = MD_TERM(num_dims,{term1_v,term1_th});

% term2 is a simple div term, defining K(th) = ZeE/m sin(th)
%
% eq1 : df/dt == div(K(th) f)     [pterm1:div(g(th)=K(th),-1, BCL=N,BCR=N)]

% dV_v = @(x,p,t,d) x.^2; %changing for MASS term
% dV_th = @(x,p,t,d) sin(x);

dV_v = @(x,p,t,d) x; %changing for MASS term
dV_th = @(x,p,t,d) sin(x);

g1 = @(x,p,t,dat) 0*x+1;
pterm1   =  MASS(g1,'','',dV_v);
term2_v = SD_TERM({pterm1});

K = @(x,p) params_si.a.Z.*params_si.e.*params_si.E.*sin(x)./params_si.a.m;
g2 = @(x,p,t,dat) -K(x,p);

pterm1 = DIV(num_dims,g2,'',-1,'N','N','','','',dV_th);
term2_th = SD_TERM({pterm1});
term2 = MD_TERM(num_dims,{term2_v,term2_th});

% dV_v = @(x,p,t,d) x.^2;
% dV_th = @(x,p,t,d) sin(x);

%%
% term3 is done using SLDG defining A(v)= sqrt(0.5*nu_par*v^2) and B(th) =
% -1
%
% eq1 :  df/dt == div(A(v) * q)        [pterm1: div (g(v)=A(v),+1, BCL=?, BCR=?)]
% eq2 :      q == A(v) * grad(f)       [pterm2: grad(g(v)=A(v),-1, BCL=D, BCR=N)]

dV_v = @(x,p,t,d) x.^2; 
dV_th = @(x,p,t,d) sin(x);

A = @(v,p) sqrt(0.5.*v.^2.*(p.nu_par(v,p.a,p.b)+ p.nu_par(v,p.a,p.b2)));
g1 = @(v,p,t,dat) A(v,p);
g2 = @(v,p,t,dat) A(v,p);

pterm1 = DIV (num_dims,g1,'',+1,'N','D','','','',dV_v);
pterm2 = GRAD(num_dims,g2,'',-1,'D','N',BCL,BCR,'',dV_v);
term3_v = SD_TERM({pterm1,pterm2});
term3 = MD_TERM(num_dims,{term3_v,[]});

% term4 is a div using B(v) = v (m_a/(m_a + m_b))nu_s
%
% eq1 :  df/dt == div(B(v) * f)       [pterm1: div(g(p)=B(v),+1, BCL=?, BCR=?]

dV_v = @(x,p,t,d) x.^2; 
dV_th = @(x,p,t,d) sin(x);

B = @(v,p) v*p.a.m.*(p.nu_s(v,p.a,p.b)./(p.a.m + p.b.m) + p.nu_s(v,p.a,p.b2)./(p.a.m + p.b2.m));
g3 = @(v,p,t,dat) -B(v,p);

pterm1 = DIV(num_dims,g3,'',-1,'D','N',BCL,'','',dV_v);

term4_v = SD_TERM({pterm1});
    %     else
term4   = MD_TERM(num_dims,{term4_v,[]});

% term5 is done combining mass and div defining C(v) = (nu_D(v)/2) 
%
% eq1 :  df/dt == div( q)   [pterm1: div (g(v)=D(v),+1, BCL=D, BCR=D]
% eq2 :      q ==  C(v) grad(f)  [pterm2: grad(g(v)=D(v),-1, BCL=N, BCR=N]

dV_v = @(x,p,t,d) x.^2; %changing for MASS term
dV_th = @(x,p,t,d) sin(x);

C = @(v,p) sqrt((p.nu_D(v,p.a,p.b) + p.nu_D(v,p.a,p.b2))/2);
g4 = @(v,p,t,dat) C(v,p);
pterm1 = MASS(g4,[],[],dV_v);
term5_v = SD_TERM({pterm1,pterm1});

D = @(v,p) v.*0 + 1;
g5 = @(v,p,t,dat) D(v,p);
pterm1 = DIV (num_dims,g5,'',+1,'D','D','','','',dV_th);
pterm2 = GRAD(num_dims,g5,'',-1,'N','N','','','',dV_th);
term5_th = SD_TERM({pterm1,pterm2});
term5   = MD_TERM(num_dims,{term5_v,term5_th});

terms = {term1b};
%% Define sources 

sources = {};

%% Define function to set time step
    function dt=set_dt(pde,CFL)    
        Lmax = pde.dimensions{1}.max;
        LevX = pde.dimensions{1}.lev;
        dt = Lmax/2^LevX*CFL;
    end


%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params_si,@set_dt,[],initial_conditions,solutions);

end
