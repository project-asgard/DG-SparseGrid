 function pde = mirror2_space_pitch_div(opts)

% Two-dimensional test of spatial and pitch-angle variation in magnetic
% mirror where f coefficients change signs throughout the domain
% 
% df/dt ==  1/(sin(th)) d/dth (sin(th) sin(th)/2B dB/ds f_a)
%           - 1/B d (B cos(th)f)/ds
% Run with
%
% asgard(@mirror2_space_pitch_div,'timestep_method','BE','case',3,'dt',1e-8,'lev',3,'deg', 3,'num_steps',20)

params = mirror_parameters();

params.B_func = @(s) s.^2;
params.dB_ds = @(s) 2.*s;
        
mag = @(s) params.dB_ds(s)./(2*params.B_func(s));

dim_th = DIMENSION(0,pi);
dV_th = @(x,p,t,d) sin(x);
dim_th.moment_dV = dV_th;

dim_s = DIMENSION(-3,3);
dV_s = @(x,p,t,d) params.B_func(x)/params.a.m;
dim_s.moment_dV = dV_s;

dimensions = {dim_th,dim_s};
num_dims = numel(dimensions);

%% Define the initial conditions

params.init_cond_z = @(z,p,t,dat) cos(z);
params.init_cond_s = @(s,p,t,dat) s;
ic1 = new_md_func(num_dims,{params.init_cond_z,params.init_cond_s});
initial_conditions = {ic1};

%% Define the boundary conditions

BCL = new_md_func(num_dims,{...
    @(x,p,t,dat) cos(x), ...
    @(s,p,t,dat) s, ...
    params.boundary_cond_t});

BCR = new_md_func(num_dims,{...
    @(x,p,t,dat) cos(x), ...
    @(s,p,t,dat) s, ...
    params.boundary_cond_t});

%% Define the terms of the PDE

%% Advection  term

%termZa is a done combining mass ad div defining C(th) = sin(th) and using
%the function mag(s) = dB/ds/(2*B). This applies to the domain where 
% mag(s) > 0 and therefore we use upwinding
%eq1: df/dt == div(mag(s) f)       [pterm1: div(g(th)=C(th),-1, BCL=N,BCR=N)]

dV_th = @(x,p,t,d) sin(x);
dV_s = @(x,p,t,d) params.B_func(x)/params.a.m;

g1 = @(x,p,t,d) mag(x).*(mag(x) > 0);
pterm1 = MASS(g1,'','',dV_s);
termZa_s = SD_TERM({pterm1});

C = @(z) sin(z);
g2 = @(z,p,t,dat) C(z);
pterm2 = DIV(num_dims,g2,'',-1,'N','D','',BCR,'',dV_th);
termZa_z = SD_TERM({pterm2});

termZa = MD_TERM(num_dims,{termZa_z,termZa_s});

%termZb is a done combining mass ad div defining C(th) = sin(th) and using
%the function mag(s) = dB/ds/(2*B). This applies to the domain where 
% mag(s) < 0 and therefore we use downwinding
%eq1: df/dt == div(mag(s) f)       [pterm1: div(g(th)=C(th),+1, BCL=N,BCR=N)]

dV_th = @(x,p,t,d) sin(x);
dV_s = @(x,p,t,d) params.B_func(x)/params.a.m;

g1 = @(x,p,t,d) mag(x).*(mag(x) < 0);
pterm1 = MASS(g1,'','',dV_s);
termZb_s = SD_TERM({pterm1});

C = @(z)  sin(z);
g2 = @(z,p,t,dat) C(z);
pterm2 = DIV(num_dims,g2,'',+1,'D','N',BCL,'','',dV_th);
termZb_z = SD_TERM({pterm2});

termZb = MD_TERM(num_dims,{termZb_z,termZb_s});

%termSa is done combining mass and div G(th) =
%cos(th) and H(s) = -1. Here we restrict to z < pi/2 and use upwinding
%eq1 : df/dt == div(G(th) f)      [pterm1: div(g(s)=H(s),-1, BCL=D,BCR=N)]

dV_th = @(x,p,t,d) sin(x);
dV_s = @(x,p,t,d) params.B_func(x)./params.a.m;

G = @(z) cos(z).*(z < pi/2);
g2 = @(z,p,t,dat) G(z);
pterm1 = MASS(g2,'','',dV_th);
termSa_z = SD_TERM({pterm1});

H = @(s) s.*0 - 1;
g3 = @(s,p,t,dat) H(s);
pterm2 = DIV(num_dims,g3,'',-1,'D','N','',BCR,'',dV_s);
termSa_s = SD_TERM({pterm2});

termSa  = MD_TERM(num_dims,{termSa_z,termSa_s});

%termSb is done combining mass and div G(th) =
%cos(th) and H(s) = -1. Here we restrict to z > pi/2 and use downwinding
%eq1 : df/dt == div(G(th) f)      [pterm1: div(g(s)=H(s),-1, BCL=D,BCR=N)]

dV_th = @(x,p,t,d) sin(x);
dV_s = @(x,p,t,d) params.B_func(x)/params.a.m;

G = @(z) cos(z).*(z > pi/2);
g2 = @(z,p,t,dat) G(z);
pterm1 = MASS(g2,'','',dV_th);
termSb_z = SD_TERM({pterm1});

H = @(s) s.*0 - 1;
g3 = @(s,p,t,dat) H(s);
pterm2 = DIV(num_dims,g3,'',+1,'N','D',BCL,'','',dV_s);
termSb_s = SD_TERM({pterm2});

termSb  = MD_TERM(num_dims,{termSb_z,termSb_s});

%% Mass term
% termM == cos(z)dB/ds/(2*B) f
% termM == r(z)*w(s)
% r(z) == g1(z) [mass, g1(z) = -cos(z),  BC N/A]
% w(s) == g2(s) f [mass, g2(s) = mag(s), BC N/A]

% dV_th = @(x,p,t,d) sin(x);
% dV_s = @(x,p,t,d) params.B_func(x);
% 
% G = @(z) cos(z);
% g2 = @(z,p,t,dat) G(z);
% pterm2 = MASS(g2,'','',dV_th);
% termB1_z = SD_TERM({pterm2});
% 
% g3 = @(s,p,t,dat) mag(s);
% pterm3 = MASS(g3,'','',dV_s);
% termB1_s = SD_TERM({pterm3});
% 
% termM = MD_TERM(num_dims,{termB1_z,termB1_s});

terms = {termZa,termZb,termSa,termSb};



%% Define the analytic solution (optional)


soln1_z = @(z,p,t,d) cos(z);
soln1_s = @(s,p,t,d) s;
solution = new_md_func(num_dims,{soln1_z,soln1_s});
solutions = {solution};


%% Define sources
%source1_v = @(x,p,t,d) p.f0_v(x);
source1_z = @(x,p,t,d) (sin(x)).^2 - (cos(x)).^2;
source1 = new_md_func(num_dims,{source1_z,[]});
sources = {source1};

%% Define function to set time step
    function dt=set_dt(pde,CFL)      
        Lmax = pde.dimensions{1}.max;
        LevX = pde.dimensions{1}.lev;
        dt = Lmax/2^LevX*CFL;
    end



%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,[],initial_conditions,solutions);
end