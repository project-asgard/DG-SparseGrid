function pde = mirror3_div(opts)
% Three-dimensional magnetic mirror from the FP paper - evolution of the ion velocity and spatial dependence
% of f in the presence of Coulomb collisions with background electrons
% using div operator
% 
% df/dt == (v^2 sin(z)^2*cos(z))/2 dB/ds/B)df/dv - vcos(z)df/ds - vcos(z)f dB/ds/B + nu_D/(2*sin(z)) d/dz ( sin(z) df/dz ) + 1/v^2 (d/dv(v^3[(m_a/(m_a + m_b))nu_s f) + 0.5*nu_par*v*d/dv(f)]))
%
% df/dt == (v^2 sin(z)^2*cos(z))/2 dB/ds/B)df/dv - vcos(z)df/ds - vcos(z)f dB/ds/B + nu_D/(2*sin(z)) d/dz ( sin(z) df/dz ) + 1/v^2 (d/dv(flux_v))
%a
% flux_v == v^3[(m_a/(m_a + m_b))nu_s f) + 0.5*nu_par*v*d/dv(f)]
%
% Run with
%
% asgard(@mirror3_div,'timestep_method','BE','case',3)

% choose a case
%
% case 1 - maxwellian offset and different Temp.
% case 2 - max. no offset and different Temp.
% case 3 - max. offset and same Temp.ans

params = mirror_parameters();

switch opts.case_
    case 1
        %initial test with no collisions and only magnetic field varying in
        %space
        params.a.Z = 1; %hydrogen particles in beam
        params.b.Z = 1; %hydrogen particles in background
        params.a.m = params.m_H; %hydrogen mass in beam
        params.b.m = params.m_H; %hydrogen mass in background
        params.a.T_eV = 10; %Target temperature in Kelvin
        params.b.E_eV = params.a.E_eV; %background energy also at 1 keV
        params.b.T_eV = params.a.T_eV;
        params.a.s0 = 0; %beam centered at s = 0
        params.a.ds0 = 0.3; %spread of beam
        params.a.dz0 = sqrt(params.a.T_eV/params.a.E_eV); %spread of pitch angle
        params.a.n = 1e19; %beam density
        params.b.n = params.a.n; %background density equal to beam density
        params.B_func = params.B_func2; %setting magnetic field to loop function
        params.dB_ds = params.dB_ds2;
        params.init_cond_v = @(v,p,t) params.a.n.*params.gauss(v,params.a.v_beam,params.a.vth);
        params.init_cond_s = @(s,p,t) params.gauss(s,params.a.s0,params.a.ds0);
    case 2
        params.init_cond_s = @(s,p,t) params.maxwell(s,params.Lx/3,2.5);
        params.boundary_cond_s = @(s,p,t) s.*0;
        params.soln_z = @(z,p,t) z.*0 + 1;
        params.soln_s = @(s,p,t) params.maxwell(s,params.coil_coords(1), 0.1) + params.maxwell(s,params.coil_coords(2),0.1);
    case 3
        params.B_func = params.B_func2; %setting magnetic field to loop function
        params.dB_ds = params.dB_ds2;
        params.init_cond_v = @(v,p,t) params.gauss(v,params.a.v_beam,params.a.vth)/params.a.n;
        params.init_cond_s = @(s,p,t) params.gauss(s,0,2);
        params.init_cond_z = @(z,p,t) z.*0 + 1;
        params.boundary_cond_s = @(s,p,t) s.*0;
        params.soln_z = @(z,p,t) z.*0 + 1;
        %params.z0 = pi/2 -1e-6;
        params.soln_s = @(s,p,t) params.maxwell(s,params.coil_coords(1), 0.1) + params.maxwell(s,params.coil_coords(2),0.1);
end

%% Define the dimensions

dim_v = DIMENSION(0,2e7);
dV_v = @(x,p,t,d) x.^2;
dim_v.moment_dV = dV_v;

dim_th = DIMENSION(0,pi);
dV_th = @(x,p,t,d) sin(x);
dim_th.moment_dV = dV_th;

dim_s = DIMENSION(-3,3);
dV_s = @(x,p,t,d) p.B_func(s); 
dim_s.moment_dV = dV_s;

dimensions = {dim_v,dim_th,dim_s};
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
pterm3 = GRAD(num_dims,g3,0,'N','N');

term1_v = SD_TERM({pterm1});
term1_z = SD_TERM({pterm2});
term1_s = SD_TERM({pterm3});
termS1  = MD_TERM(num_dims,{term1_v,term1_z,term1_s});

%% Mass term
% termS2 == -vcos(z)dB/ds/B f
% termS1 == q(v)*r(z)*w(s)
% q(v) == g1(v)  [mass, g1(p) = v,  BC N/A]
% r(z) == g2(z) [mass, g2(z) = cos(z),  BC N/A]
% w(s) == g3(s) f [mass, g3(s) = -dB/ds/B, BC N/A]

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

% term V1 == (sin(z)^2*cos(z)/2 *dB/ds/B)v^2 df/dv
% term V1 == w(v)q(z)r(s)
% w(v) = g1(v) m(v) [mass, g1(v) = v.^2, BC = N/A]
% m(v) = d/dv (g2(v) f) [grad, g2(v) = 1, BCL=N, BCR=D ]
% q(z) = g3(z) [mass, g3(z) = 0.5*sin(z)^2 cos(z), BC = N/A]
% r(s) = g4(s) [mass, g5(s) = dB/ds / B, BC = N/A]
g1 = @(v,p,t,dat) v.^2;
g2 = @(v,p,t,dat) v.*0 + 1;
pterm1 = MASS(g1);
pterm2 = GRAD(num_dims,g2,0,'N','N');%, BCL, BCR);
termV_v = SD_TERM({pterm1,pterm2});

g3 = @(z,p,t,dat) 0.5*sin(z).^2.*cos(z);
pterm1 = MASS(g3);
termV_z = SD_TERM({pterm1});

g4 = @(s,p,t,dat) p.dB_ds(s)./p.B_func(s);
pterm1 = MASS(g4);
termV_s = SD_TERM({pterm1});
termV1 = MD_TERM(num_dims,{termV_v,termV_z,termV_s});
%termV1 = MD_TERM(num_dims,{termV_v,[],termV_s});

terms = {termV1,termS1,termS2,termC1,termC2,termC3};

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
