function pde = continuity6
% 3D test case using continuity equation, i.e.,
% df/dt + v.grad_x(f) + a.grad_v(f)==0 where v={1,1,3}, a={4,3,2}

%% Setup the dimensions
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

vx=1;
vy=1;
vz=3;
ax=4;
ay=3;
az=2;

%
% Here we setup a 6D problem (x,y,z,vx,vy,vz)

dim_x.name = 'x';
dim_x.BCL = 0; % periodic
dim_x.BCR = 0;
dim_x.domainMin = -1;
dim_x.domainMax = +1;
dim_x.lev = 2;
dim_x.deg = 2;
dim_x.FMWT = []; % Gets filled in later
dim_x.init_cond_fn = @(x,p) x.*0;

dim_y.name = 'y';
dim_y.BCL = 0; % periodic
dim_y.BCR = 0;
dim_y.domainMin = -2;
dim_y.domainMax = +2;
dim_y.lev = 2;
dim_y.deg = 2;
dim_y.FMWT = []; % Gets filled in later
dim_y.init_cond_fn = @(y,p) y.*0;

dim_z.name = 'z';
dim_z.BCL = 0; % periodic
dim_z.BCR = 0;
dim_z.domainMin = -3;
dim_z.domainMax = +3;
dim_z.lev = 2;
dim_z.deg = 2;
dim_z.FMWT = []; % Gets filled in later
dim_z.init_cond_fn = @(z,p) z.*0;

dim_vx.name = 'vx';
dim_vx.BCL = 0; % periodic
dim_vx.BCR = 0;
dim_vx.domainMin = -10;
dim_vx.domainMax = +10;
dim_vx.lev = 2;
dim_vx.deg = 2;
dim_vx.FMWT = []; % Gets filled in later
dim_vx.init_cond_fn = @(x,p) x.*0;

dim_vy.name = 'vy';
dim_vy.BCL = 0; % periodic
dim_vy.BCR = 0;
dim_vy.domainMin = -20;
dim_vy.domainMax = +20;
dim_vy.lev = 2;
dim_vy.deg = 2;
dim_vy.FMWT = []; % Gets filled in later
dim_vy.init_cond_fn = @(y,p) y.*0;

dim_vz.name = 'vz';
dim_vz.BCL = 0; % periodic
dim_vz.BCR = 0;
dim_vz.domainMin = -30;
dim_vz.domainMax = +30;
dim_vz.lev = 2;
dim_vz.deg = 2;
dim_vz.FMWT = []; % Gets filled in later
dim_vz.init_cond_fn = @(z,p) z.*0;

%%
% Add dimensions to the pde object
% Note that the order of the dimensions must be consistent with this across
% the remainder of this PDE.

pde.dimensions = {dim_x,dim_y,dim_z,dim_vx,dim_vy,dim_vz};

%% Setup the terms of the PDE
%
% Here we have 6 terms, having only nDims=6 (x,y,z,vx,vy,vz) operators.

%%
% Setup the v_x.d_dx (v . GradX . MassY . MassZ) term

term2_x.type = 1; % grad (see coeff_matrix.m for available types)
term2_x.G = @(x,t,dat) x*0+vx; 
term2_x.TD = 0; 
term2_x.dat = []; 
term2_x.LF = 0; 
term2_x.name = 'v_x.d_dx';

term2 = term_fill({term2_x,[],[],[],[],[]});

%%
% Setup the v_y.d_dy (v . MassX . GradY . MassZ) term

term3_y.type = 1; % grad (see coeff_matrix.m for available types)
term3_y.G = @(y,t,dat) y*0+vy; 
term3_y.TD = 0; 
term3_y.dat = []; 
term3_y.LF = 0; 
term3_y.name = 'v_y.d_dy';

term3 = term_fill({[],term3_y,[],[],[],[]});

%%
% Setup the v_z.d_dz (v . MassX . MassY . GradZ) term

term4_z.type = 1; % grad (see coeff_matrix.m for available types)
term4_z.G = @(z,t,dat) z*0+vz; 
term4_z.TD = 0; 
term4_z.dat = []; 
term4_z.LF = 0; 
term4_z.name = 'v_z.d_dz';

term4 = term_fill({[],[],term4_z,[],[],[]});

%%
% Setup the a_x.d_dvx (a . GradVX . MassVY . MassVZ) term

term5_vx.type = 1; % grad (see coeff_matrix.m for available types)
term5_vx.G = @(x,t,dat) x*0+ax; 
term5_vx.TD = 0; 
term5_vx.dat = []; 
term5_vx.LF = 0; 
term5_vx.name = 'a_x.d_dvx';

term5 = term_fill({[],[],[],term5_vx,[],[]});

%%
% Setup the a_y.d_dvy (a . MassVX . GradVY . MassVZ) term

term6_vy.type = 1; % grad (see coeff_matrix.m for available types)
term6_vy.G = @(y,t,dat) y*0+ay; 
term6_vy.TD = 0; 
term6_vy.dat = []; 
term6_vy.LF = 0; 
term6_vy.name = 'a_y.d_dvy';

term6 = term_fill({[],[],[],[],term6_vy,[]});

%%
% Setup the a_z.d_dvz (a . MassVX . MassVY . GradVZ) term

term7_vz.type = 1; % grad (see coeff_matrix.m for available types)
term7_vz.G = @(z,t,dat) z*0+az; 
term7_vz.TD = 0; 
term7_vz.dat = []; 
term7_vz.LF = 0; 
term7_vz.name = 'a_z.d_dvz';

term7 = term_fill({[],[],[],[],[],term7_vz});

%%
% Add terms to the pde object

pde.terms = {term2,term3,term4,term5,term6,term7};

%% Construct some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

pde.params = params;

%% Add an arbitrary number of sources to the RHS of the PDE
% Each source term must have nDims + 1 (which here is 2+1 / (x,v) + time) functions describing the 
% variation of each source term with each dimension and time.
% Here we define 3 source terms.

%%
% Source term 0
source0 = { ...
    @(x,p) cos(pi*x),         ... % s1x
    @(y,p) sin(2*pi*y),       ... % s1y
    @(z,p) cos(2*pi*z/3),     ... % s1z
    @(vx,p) cos(pi*vx/5),     ... % s1vx
    @(vy,p) sin(3*pi*vy/20),  ... % s1vy
    @(vz,p) cos(pi*vz/30),    ... % s1vz
    @(t)   2*cos(2*t)         ... % s1t
    };

%%
% Source term 1
source1 = { ...
    @(x,p) cos(pi*x),      ... % s2x
    @(y,p) cos(2*pi*y),    ... % s2y
    @(z,p) cos(2*pi*z/3),  ... % s2z
    @(vx,p) cos(pi*vx/5),    ... % s2vx
    @(vy,p) sin(3*pi*vy/20),  ... % s2vy
    @(vz,p) cos(pi*vz/30),... % s2vz
    @(t)   2*pi*sin(2*t)   ... % s2t
    };

%%
% Source term 2
source2 = { ...
    @(x,p) sin(pi*x),     ... % s3x
    @(y,p) sin(2*pi*y),   ... % s3y
    @(z,p) cos(2*pi*z/3), ... % s3z
    @(vx,p) cos(pi*vx/5),     ... % s3vx
    @(vy,p) sin(3*pi*vy/20),   ... % s3vy
    @(vz,p) cos(pi*vz/30), ... % s3vz
    @(t)   -pi*sin(2*t)   ... % s3t
    };

%%
% Source term 3
source3 = { ...
    @(x,p) cos(pi*x),       ... % s4x
    @(y,p) sin(2*pi*y),     ... % s4y
    @(z,p) sin(2*pi*z/3),   ... % s4z
    @(vx,p) cos(pi*vx/5),       ... % s4vx
    @(vy,p) sin(3*pi*vy/20),     ... % s4vy
    @(vz,p) cos(pi*vz/30),   ... % s4vz
    @(t)   -2*pi*sin(2*t) ... % s4t
    };

%%
% Source term 4
source4 = { ...
    @(x,p) cos(pi*x),       ... % s4x
    @(y,p) sin(2*pi*y),     ... % s4y
    @(z,p) cos(2*pi*z/3),   ... % s4z
    @(vx,p) cos(pi*vx/5),       ... % s4vx
    @(vy,p) cos(3*pi*vy/20),     ... % s4vy
    @(vz,p) cos(pi*vz/30),   ... % s4vz
    @(t)   9/20*pi*sin(2*t) ... % s4t
    };

%%
% Source term 5
source5 = { ...
    @(x,p) cos(pi*x),       ... % s4x
    @(y,p) sin(2*pi*y),     ... % s4y
    @(z,p) cos(2*pi*z/3),   ... % s4z
    @(vx,p) cos(pi*vx/5),       ... % s4vx
    @(vy,p) sin(3*pi*vy/20),     ... % s4vy
    @(vz,p) cos(pi*vz/30),   ... % s4vz
    @(t)   -4/5*pi*sin(2*t) ... % s4t
    };

%%
% Source term 6
source6 = { ...
    @(x,p) cos(pi*x),       ... % s4x
    @(y,p) sin(2*pi*y),     ... % s4y
    @(z,p) cos(2*pi*z/3),   ... % s4z
    @(vx,p) cos(pi*vx/5),       ... % s4vx
    @(vy,p) sin(3*pi*vy/20),     ... % s4vy
    @(vz,p) sin(pi*vz/30),   ... % s4vz
    @(t)   -1/15*pi*sin(2*t) ... % s4t
    };

%%
% Add sources to the pde data structure
pde.sources = {source0,source1,source2,source3,source4,source5,source6};

%% Define the analytic solution (optional).
% This requires nDims+time function handles.

pde.analytic_solutions_1D = { ...
    @(x,p) cos(pi*x),     ... % a_x
    @(y,p) sin(2*pi*y),   ... % a_y
    @(z,p) cos(2*pi*z/3), ... % a_z
    @(t)   sin(2*t)       ... % a_t
    };

%% Other workflow options that should perhpas not be in the PDE?

pde.set_dt = @set_dt; % Function which accepts the pde (after being updated with CMD args).
pde.solvePoisson = 0; % Controls the "workflow" ... something we still don't know how to do generally.
pde.applySpecifiedE = 0; % Controls the "workflow" ... something we still don't know how to do generally.
pde.implicit = 0; % Can likely be removed and be a runtime argument.
pde.checkAnalytic = 1; % Will only work if an analytic solution is provided within the PDE.

end

%% Define the various input functions specified above.

% function f=f0_x(x,p); f=x.*0; end
% function f=f0_y(y,p); f=y.*0; end
% function f=f0_z(z,p); f=z.*0; end

%%
% Function to set time step
function dt=set_dt(pde)

Lmax = pde.dimensions{1}.domainMax;
LevX = pde.dimensions{1}.lev;
CFL = pde.CFL;

dt = Lmax/2^LevX*CFL;
end
