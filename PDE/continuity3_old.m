function pde = continuity3_old
% 3D test case using continuity equation, i.e.,
%
% df/dt + v.grad(f)==0 where v={1,1,1}, so 
%
% df/dt = -df/dx -df/dy - df/dz
%
% Run with 
%
% explicit
% asgard(continuity3);
%
% implicit
% asgard(continuity3,'timestep_method','CN');

%% Setup the dimensions

%
% Here we setup a 3D problem (x,y,z)

dim_x.name = 'x';
dim_x.BCL = 'P'; % periodic
dim_x.BCR = 'P';
dim_x.domainMin = -1;
dim_x.domainMax = +1;
dim_x.init_cond_fn = @(x,p,t) x.*0;
dim_x.lev = 5;

dim_y.name = 'y';
dim_y.BCL = 'P'; % periodic
dim_y.BCR = 'P';
dim_y.domainMin = -2;
dim_y.domainMax = +2;
dim_y.init_cond_fn = @(y,p,t) y.*0;
dim_y.lev = 4;

dim_z.name = 'z';
dim_z.BCL = 'P'; % periodic
dim_z.BCR = 'P';
dim_z.domainMin = -3;
dim_z.domainMax = +3;
dim_z.init_cond_fn = @(z,p,t) z.*0;
dim_z.lev = 3;

%%
% Add dimensions to the pde object
% Note that the order of the dimensions must be consistent with this across
% the remainder of this PDE.

pde.dimensions = {dim_x,dim_y,dim_z};

%% Setup the terms of the PDE
%
% Here we have 3 terms, having only nDims=3 (x,y,z) operators.

%%
% -df/dx

term2_x.type = 'grad';
term2_x.G = @(x,p,t,dat) x*0+1; % G function for use in coeff_matrix construction.
term2_x.LF = 0; % central flux

term2 = term_fill({term2_x,[],[]});

%%
% -df/fy

term3_y.type = 'grad'; % grad (see coeff_matrix.m for available types)
term3_y.G = @(y,p,t,dat) y*0+1; % G function for use in coeff_matrix construction.
term3_y.LF = 0; % central flux

term3 = term_fill({[],term3_y,[]});

%%
% -df/dz

term4_z.type = 'grad'; % grad (see coeff_matrix.m for available types)
term4_z.G = @(z,p,t,dat) z*0+1; % G function for use in coeff_matrix construction.
term4_z.LF = 0; % central flux

term4 = term_fill({[],[],term4_z});

%%
% Add terms to the pde object

pde.terms = {term2,term3,term4};

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
% Source term 1
source1 = { ...
    @(x,p,t) cos(pi*x),     ... % s1x
    @(y,p,t) sin(2*pi*y),   ... % s1y
    @(z,p,t) cos(2*pi*z/3), ... % s1z
    @(t)   2*cos(2*t)     ... % s1t
    };

%%
% Source term 2
source2 = { ...
    @(x,p,t) cos(pi*x),     ... % s2x
    @(y,p,t) cos(2*pi*y),   ... % s2y
    @(z,p,t) cos(2*pi*z/3), ... % s2z
    @(t)   2*pi*sin(2*t)  ... % s2t
    };

%%
% Source term 3
source3 = { ...
    @(x,p,t) sin(pi*x),     ... % s3x
    @(y,p,t) sin(2*pi*y),   ... % s3y
    @(z,p,t) cos(2*pi*z/3), ... % s3z
    @(t)   -pi*sin(2*t)   ... % s3t
    };

%%
% Source term 4
source4 = { ...
    @(x,p,t) cos(pi*x),       ... % s4x
    @(y,p,t) sin(2*pi*y),     ... % s4y
    @(z,p,t) sin(2*pi*z/3),   ... % s4z
    @(t)   -2/3*pi*sin(2*t) ... % s4t
    };

%%
% Add sources to the pde data structure
pde.sources = {source1,source2,source3,source4};

%% Define the analytic solution (optional).
% This requires nDims+time function handles.

pde.analytic_solutions_1D = { ...
    @(x,p,t) cos(pi*x),     ... % a_x
    @(y,p,t) sin(2*pi*y),   ... % a_y
    @(z,p,t) cos(2*pi*z/3), ... % a_z
    @(t)   sin(2*t)       ... % a_t
    };

%%
% Function to set time step

    function dt=set_dt(pde,CFL)    
        Lmax = pde.dimensions{1}.domainMax;
        Lmin = pde.dimensions{1}.domainMin;
        LevX = pde.dimensions{1}.lev;
        dt = (Lmax-Lmin)/2^LevX*CFL;
    end

pde.set_dt = @set_dt;

end


