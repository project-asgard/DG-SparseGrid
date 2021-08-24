function pde = mirror2_advection_Efield_div(opts)

params_si = mirror_parameters();


switch opts.case_
    case 1 
        params_si.a.T_eV = 0.05*params_si.b.T_eV; %Target temperature in Kelvin\
        offset = 10^6; %case with offset and change in Temperature
    case 2 
        params_si.a.T_eV = 0.05*params_si.b.T_eV;
        offset = 0; %case with no offset but change in Temperature
    case 3 
<<<<<<< HEAD
        n_cgs = 5e13; %equilibrium density in cm.^-3
=======
        n_cgs = 8e14; %equilibrium density in cm.^-3
>>>>>>> Changed sign of G function in term1b inside mirror2_velocity_div.m. Code appears to be relatively stable up to dt=1e-5. Also created PDE/mirror2_advection_Efield_div.m to test out advection terms from electric field acceleration.
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
<<<<<<< HEAD
        params_si.a.T_eV = 1000;%2/3*params_si.a.E_eV;
=======
        params_si.a.T_eV = 2/3*params_si.a.E_eV;
>>>>>>> Changed sign of G function in term1b inside mirror2_velocity_div.m. Code appears to be relatively stable up to dt=1e-5. Also created PDE/mirror2_advection_Efield_div.m to test out advection terms from electric field acceleration.
        params_si.b.T_eV = params_si.a.T_eV;
        params_si.b2.T_eV = params_si.a.T_eV;
        params_si.a.vth = params_si.v_th(params_si.a.T_eV,params_si.a.m);
        params_si.b.vth = params_si.v_th(params_si.b.T_eV,params_si.b.m);
        params_si.b2.v.th = params_si.v_th(params_si.b2.T_eV,params_si.b2.m);
        params_si.ln_delt = 15;
        E_dreicer_si = params_si.a.n.*params_si.e^3*params_si.ln_delt/(4*pi*params_si.eps0^2*params_si.a.m ... 
            *params_si.a.vth^2);
<<<<<<< HEAD
        params_si.E = 10^-3*E_dreicer_si;
        %vel_norm = @(v,vth) v./vth; %normalized velocity to thermal velocity
        params_si.maxwell = @(v,offset,vth) params_si.a.n/(pi.^(3/2)*vth^3).*exp(-((v-offset)/vth).^2);
        params_si.init_cond_v = @(v,p,t) params_si.maxwell(v,0,params_si.a.vth);
        params_si.soln_v = @(v,p,t) solution_v(v,p,t);
=======
        params_si.E = 10^-6*E_dreicer_si;
        %vel_norm = @(v,vth) v./vth; %normalized velocity to thermal velocity
        params_si.maxwell = @(v,offset,vth) params_si.a.n/(pi.^(3/2)*vth^3).*exp(-((v-offset)/vth).^2);
        params_si.init_cond_v = @(v,p,t) params_si.maxwell(v,0,params_si.a.vth);
>>>>>>> Changed sign of G function in term1b inside mirror2_velocity_div.m. Code appears to be relatively stable up to dt=1e-5. Also created PDE/mirror2_advection_Efield_div.m to test out advection terms from electric field acceleration.
        %params_cgs.nu_ab0  = @(a,b) b.n * params_cgs.e^4 * a.Z^2 * b.Z^2 * params_cgs.ln_delt / (pi^3/2.*a.m^2*b.vth^3); %scaling coefficient
        %params.eps0 = 1/(4*pi);
end

<<<<<<< HEAD
    function ret = solution_v(v,p,t)
        ret = params_si.a.n/(pi^3/2.*params_si.v_th(params_si.b.T_eV,params_si.a.m).^3).*...
            exp(-(v./params_si.v_th(params_si.b.T_eV,params_si.a.m)).^2);
        if isfield(p,'norm_fac')
            ret = p.norm_fac .* ret;
        end
    end

dim_v = DIMENSION(0,15*params_si.a.vth);
=======
dim_v = DIMENSION(0,10*params_si.a.vth);
>>>>>>> Changed sign of G function in term1b inside mirror2_velocity_div.m. Code appears to be relatively stable up to dt=1e-5. Also created PDE/mirror2_advection_Efield_div.m to test out advection terms from electric field acceleration.
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

<<<<<<< HEAD
% termE1a is done combining mass and div defining F(th) = (ZeE/m cos(th)) and
=======
% term1 is done combining mass and div defining F(th) = (ZeE/m cos(th)) and
>>>>>>> Changed sign of G function in term1b inside mirror2_velocity_div.m. Code appears to be relatively stable up to dt=1e-5. Also created PDE/mirror2_advection_Efield_div.m to test out advection terms from electric field acceleration.
% G(v) = 1
%
% eq1 : df/dt == div(F(th) f)      [pterm1: div(g(v)=G(v),-1, BCL=D,BCR=N)]

<<<<<<< HEAD
 dV_v = @(x,p,t,d) x.^2; 
 dV_th = @(x,p,t,d) sin(x);

F = @(x,p) cos(x);
g1 = @(x,p,t,dat) F(x,p).*(x>pi/2);
pterm1 = MASS(g1,'','',dV_th);
termE1_th = SD_TERM({pterm1});

G = @(v,p) v.*0 - params_si.a.Z.*params_si.e.*params_si.E./params_si.a.m;
g2 = @(v,p,t,dat) G(v,p);
pterm1 = DIV(num_dims,g2,'',+1,'D','N',BCL,'','',dV_v);
termE1_v = SD_TERM({pterm1});
termE1a = MD_TERM(num_dims,{termE1_v,termE1_th});

%termE1b is the same form as term1 but accounting for the flow in the
=======
% dV_v = @(x,p,t,d) x; %changing for MASS term
% dV_th = @(x,p,t,d) sin(x);

F = @(x,p) cos(x);
g1 = @(x,p,t,dat) F(x,p).*(x>pi/2);
pterm1 = MASS(g1,[],[],dV_th);
term1_th = SD_TERM({pterm1});

G = @(v,p) v.*0 - 1;
g2 = @(v,p,t,dat) G(v,p);
pterm1 = DIV(num_dims,g2,'',+1,'D','N',BCL,'','',dV_v);
term1_v = SD_TERM({pterm1});
term1a = MD_TERM(num_dims,{term1_v,term1_th});

%term1b is the same form as term1 but accounting for the flow in the
>>>>>>> Changed sign of G function in term1b inside mirror2_velocity_div.m. Code appears to be relatively stable up to dt=1e-5. Also created PDE/mirror2_advection_Efield_div.m to test out advection terms from electric field acceleration.
%opposite direction

F = @(x,p) cos(x);
g1 = @(x,p,t,dat) F(x,p).*(x<pi/2);
<<<<<<< HEAD
pterm1 = MASS(g1,'','',dV_th);
termE1_th = SD_TERM({pterm1});

G = @(v,p) v.*0 - params_si.a.Z.*params_si.e.*params_si.E./params_si.a.m;
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

g1 = @(x,p,t,dat) 0*x+1;
pterm1   =  MASS(g1,'','',dV_v);
termE2_v = SD_TERM({pterm1});

K = @(x,p) params_si.a.Z.*params_si.e.*params_si.E.*sin(x)./params_si.a.m;
g2 = @(x,p,t,dat) K(x,p);

pterm1 = DIV(num_dims,g2,'',-1,'N','N',BCL,BCR,'',dV_th);
termE2_th = SD_TERM({pterm1});
termE2 = MD_TERM(num_dims,{termE2_v,termE2_th});

terms = {termE1a,termE1b,termE2};
=======
pterm1 = MASS(g1,[],[],dV_th);
term1_th = SD_TERM({pterm1});

G = @(v,p) v.*0 + 5e9;
g2 = @(v,p,t,dat) G(v,p);
pterm1 = DIV(num_dims,g2,'',-1,'N','D','',BCR,'',dV_v);
term1_v = SD_TERM({pterm1});
term1b = MD_TERM(num_dims,{term1_v,term1_th});

terms = {term1b};
>>>>>>> Changed sign of G function in term1b inside mirror2_velocity_div.m. Code appears to be relatively stable up to dt=1e-5. Also created PDE/mirror2_advection_Efield_div.m to test out advection terms from electric field acceleration.

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