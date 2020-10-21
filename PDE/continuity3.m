function pde = continuity3(opts)
% 3D test case using continuity equation, i.e.,
%
% df/dt + v.grad(f)==0 where v={1,1,1}, so 
%
% df/dt = -df/dx -df/dy - df/dz
%
% Run with 
%
% explicit
% asgard(@continuity3,'lev',[5,4,3],'calculate_mass',false,'quiet',true);
%
% implicit
% asgard(@continuity3,'timestep_method','CN','calculate_mass',false,'quiet',true);
%
% NOTES
% DLG - should be able to turn the mass calculation back on after we merge
% in the fix from mirror_3D (also will fix the plotting issue)

%% Define the dimensions

% Here we setup a 3D problem (x,y,z)

dim_x = DIMENSION(-1,+1);
dim_x.name = 'x';
dim_x.init_cond_fn = @(x,p,t) x.*0;
dim_x.lev = 5;

dim_y = DIMENSION(-2,+2);
dim_y.name = 'y';
dim_y.init_cond_fn = @(y,p,t) y.*0;
dim_y.lev = 4;

dim_z = DIMENSION(-3,+3);
dim_z.name = 'z';
dim_z.init_cond_fn = @(z,p,t) z.*0;
dim_z.lev = 3;

dimensions = {dim_x,dim_y,dim_z};
num_dims = numel(dimensions);

%% Define the terms of the PDE
%
% Here we have 3 terms, having only nDims=3 (x,y,z) operators.

%%
% -df/dx

g1 = @(x,p,t,dat) x.*0-1;
pterm1  = GRAD(num_dims,g1,0,'P','P');
term1_x = TERM_1D({pterm1});
term1   = TERM_ND(num_dims,{term1_x,[],[]});

%%
% -df/fy

g1 = @(y,p,t,dat) y.*0-1;
pterm1  = GRAD(num_dims,g1,0,'P','P');
term2_y = TERM_1D({pterm1});
term2   = TERM_ND(num_dims,{[],term2_y,[]});

%%
% -df/dz

g1 = @(z,p,t,dat) z*0-1;
pterm1  = GRAD(num_dims,g1,0,'P','P');
term3_z = TERM_1D({pterm1});
term3   = TERM_ND(num_dims,{[],[],term3_z});

terms = {term1,term2,term3};

%% Define some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;


%% Define sources

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

sources = {source1,source2,source3,source4};

%% Define the analytic solution (optional).
% This requires nDims+time function handles.

analytic_solutions_1D = { ...
    @(x,p,t) cos(pi*x),     ... % a_x
    @(y,p,t) sin(2*pi*y),   ... % a_y
    @(z,p,t) cos(2*pi*z/3), ... % a_z
    @(t)   sin(2*t)       ... % a_t
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


