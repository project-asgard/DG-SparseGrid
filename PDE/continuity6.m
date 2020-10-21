function pde = continuity6(opts)
% 3D test case using continuity equation, i.e.,
%
% df/dt + b.grad_x(f) + a.grad_v(f)==0 where b={1,1,3}, a={4,3,2}
%
% df/dt == -bx*df/dx - by*df/dy - bz*df/dz - ax*df/dvx - ay*df/dvy - az*df/dvz
%
% Run with 
%
% explicit
% asgard(@continuity6,'lev',2,'deg',3,'num_steps',1,'calculate_mass',false);
%
% NOTES
% DLG - seems to be failing on my laptop due to MAXSIZE restrictions?

bx=1;
by=1;
bz=3;
ax=4;
ay=3;
az=2;

%% Define the dimensions

%%
% Here we setup a 6D problem (x,y,z,vx,vy,vz)

dim_x = DIMENSION(-1,+1);
dim_x.init_cond_fn = @(x,p,t) x.*0;

dim_y = DIMENSION(-2,+2);
dim_y.init_cond_fn = @(y,p,t) y.*0;

dim_z = DIMENSION(-3,+3);
dim_z.init_cond_fn = @(z,p,t) z.*0;

dim_vx = DIMENSION(-10,+10);
dim_vx.init_cond_fn = @(x,p,t) x.*0;

dim_vy = DIMENSION(-20,+20);
dim_vy.init_cond_fn = @(y,p,t) y.*0;

dim_vz = DIMENSION(-30,+30);
dim_vz.init_cond_fn = @(z,p,t) z.*0;

dimensions = {dim_x,dim_y,dim_z,dim_vx,dim_vy,dim_vz};
num_dims = numel(dimensions);


%% Define the terms of the PDE
%
% Here we have 6 terms, having only nDims=6 (x,y,z,vx,vy,vz) operators.

%%
% -bx*df/dx

% term2_x.type = 'grad';
% term2_x.G = @(x,p,t,dat) x*0-bx; 
% term2 = term_fill({term2_x,[],[],[],[],[]});

g1 = @(x,p,t,dat) x*0-bx; 
pterm1  = GRAD(num_dims,g1,0,'P','P');
term1_x = TERM_1D({pterm1});
term1   = TERM_ND(num_dims,{term1_x,[],[],[],[],[]});


%%
% -by*df/dy

% term3_y.type = 'grad'; 
% term3_y.G = @(y,p,t,dat) y*0-by; 
% term3 = term_fill({[],term3_y,[],[],[],[]});

g1 = @(y,p,t,dat) y*0-by; 
pterm1  = GRAD(num_dims,g1,0,'P','P');
term2_y = TERM_1D({pterm1});
term2   = TERM_ND(num_dims,{[],term2_y,[],[],[],[]});

%%
% -bz*df/dz

% term4_z.type = 'grad';
% term4_z.G = @(z,p,t,dat) z*0-bz; 
% term4 = term_fill({[],[],term4_z,[],[],[]});

g1 = @(z,p,t,dat) z*0-bz;
pterm1  = GRAD(num_dims,g1,0,'P','P');
term3_z = TERM_1D({pterm1});
term3   = TERM_ND(num_dims,{[],[],term3_z,[],[],[]});

%%
% -ax*df/dvx

% term5_vx.type = 'grad'; 
% term5_vx.G = @(x,p,t,dat) x*0-ax; 
% term5 = term_fill({[],[],[],term5_vx,[],[]});

g1 = @(x,p,t,dat) x*0-ax;
pterm1   = GRAD(num_dims,g1,0,'P','P');
term4_vx = TERM_1D({pterm1});
term4    = TERM_ND(num_dims,{[],[],[],term4_vx,[],[]});

%%
% -ay*df/dvy

% term6_vy.type = 'grad'; 
% term6_vy.G = @(y,p,t,dat) y*0-ay; 
% term6 = term_fill({[],[],[],[],term6_vy,[]});

g1 = @(y,p,t,dat) y*0-ay;
pterm1   = GRAD(num_dims,g1,0,'P','P');
term5_vy = TERM_1D({pterm1});
term5    = TERM_ND(num_dims,{[],[],[],[],term5_vy,[]});

%%
% -az*df/dvz

% term7_vz.type = 'grad'; 
% term7_vz.G = @(z,p,t,dat) z*0-az; 
% term7 = term_fill({[],[],[],[],[],term7_vz});

g1 = @(z,p,t,dat) z*0-az;
pterm1   = GRAD(num_dims,g1,0,'P','P');
term6_vz = TERM_1D({pterm1});
term6    = TERM_ND(num_dims,{[],[],[],[],[],term6_vz});

terms = {term1,term2,term3,term4,term5,term6};

%% Define some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;


%% Define sources

targ = 2;
xarg = pi;
yarg = pi/2;
zarg = pi/3;
vxarg = pi/10;
vyarg = pi/20;
vzarg = pi/30;

%%
% Source term 0
source0 = { ...
    @(x,p,t)  cos(xarg*x), ... 
    @(y,p,t)  sin(yarg*y), ... 
    @(z,p,t)  cos(zarg*z), ... 
    @(vx,p,t) cos(vxarg*vx), ...
    @(vy,p,t) sin(vyarg*vy), ...
    @(vz,p,t) cos(vzarg*vz), ...
    @(t)  2*cos(targ*t)    ...
    };

%%
% Source term 1
source1 = { ...
    @(x,p,t)  cos(xarg*x), ... 
    @(y,p,t)  cos(yarg*y), ... 
    @(z,p,t)  cos(zarg*z), ... 
    @(vx,p,t) cos(vxarg*vx), ...
    @(vy,p,t) sin(vyarg*vy), ...
    @(vz,p,t) cos(vzarg*vz), ...
    @(t)  1/2*pi*sin(targ*t)    ...
    };

%%
% Source term 2
source2 = { ...
    @(x,p,t)  sin(xarg*x), ... 
    @(y,p,t)  sin(yarg*y), ... 
    @(z,p,t)  cos(zarg*z), ... 
    @(vx,p,t) cos(vxarg*vx), ...
    @(vy,p,t) sin(vyarg*vy), ...
    @(vz,p,t) cos(vzarg*vz), ...
    @(t)  -pi*sin(targ*t)    ...
    };

%%
% Source term 3
source3 = { ...
    @(x,p,t)  cos(xarg*x), ... 
    @(y,p,t)  sin(yarg*y), ... 
    @(z,p,t)  sin(zarg*z), ... 
    @(vx,p,t) cos(vxarg*vx), ...
    @(vy,p,t) sin(vyarg*vy), ...
    @(vz,p,t) cos(vzarg*vz), ...
    @(t)  -pi*sin(targ*t)    ...
    };

%%
% Source term 4
source4 = { ...
    @(x,p,t)  cos(xarg*x), ... 
    @(y,p,t)  sin(yarg*y), ... 
    @(z,p,t)  cos(zarg*z), ... 
    @(vx,p,t) cos(vxarg*vx), ...
    @(vy,p,t) cos(vyarg*vy), ...
    @(vz,p,t) cos(vzarg*vz), ...
    @(t)  3/20*pi*sin(targ*t)    ...
    };

%%
% Source term 5
source5 = { ...
    @(x,p,t)  cos(xarg*x), ... 
    @(y,p,t)  sin(yarg*y), ... 
    @(z,p,t)  cos(zarg*z), ... 
    @(vx,p,t) sin(vxarg*vx), ...
    @(vy,p,t) sin(vyarg*vy), ...
    @(vz,p,t) cos(vzarg*vz), ...
    @(t)  -2/5*pi*sin(targ*t)    ...
    };

%%
% Source term 6
source6 = { ...
    @(x,p,t)  cos(xarg*x), ... 
    @(y,p,t)  sin(yarg*y), ... 
    @(z,p,t)  cos(zarg*z), ... 
    @(vx,p,t) cos(vxarg*vx), ...
    @(vy,p,t) sin(vyarg*vy), ...
    @(vz,p,t) sin(vzarg*vz), ...
    @(t)  -1/15*pi*sin(targ*t)    ...
    };

sources = {source0,source1,source2,source3,source4,source5,source6};

%% Define the analytic solution (optional).
% This requires nDims+time function handles.

analytic_solutions_1D = { ...
    @(x,p,t) cos(xarg*x), ...
    @(y,p,t) sin(yarg*y), ... 
    @(z,p,t) cos(zarg*z), ... 
    @(vx,p,t) cos(vxarg*vx), ... 
    @(vy,p,t) sin(vyarg*vy), ... 
    @(vz,p,t) cos(vzarg*vz), ...
    @(t)   sin(targ*t) 
    };

%% Define function to set time step
    function dt=set_dt(pde,CFL)     
        Lmax = pde.dimensions{1}.max;
        Lmin = pde.dimensions{1}.min;
        LevX = pde.dimensions{1}.lev;        
        dt = (Lmax-Lmin)/2^LevX*CFL;
    end

%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,analytic_solutions_1D);

end


