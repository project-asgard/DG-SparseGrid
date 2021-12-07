function pde = continuity2(opts)
% 2D test case using continuity equation, i.e.,
%
% df/dt == -df/dx - df/dy
%
% Run with
%
% explicit
% asgard(@continuity2_div,'lev',3,'deg',3,'CFL',0.1)
%
% implicit
% asgard(@continuity2_div,'timestep_method','CN')
%
% with adaptivity
% asgard(@continuity2_div,'timestep_method','CN','adapt',true)

% pde.CFL = 0.1;

soln_x = @(x,p,t)  cos(pi*x);
soln_y = @(y,p,t)  sin(2*pi*y);
soln_t = @(t)      sin(2*t);

%% Define the dimensions
%
% Here we setup a 2D problem (x,y)

dim_x = DIMENSION(-1,+1);
dim_x.moment_dV = @(x,p,t,dat) 0*x+1;
dim_y = DIMENSION(-2,+2);
dim_y.moment_dV = @(y,p,t,dat) 0*y+1;
dimensions = {dim_x,dim_y};
num_dims = numel(dimensions);

%% Initial conditions
ic1 = new_md_func(num_dims,{soln_x,soln_y,soln_t});
initial_conditions = {ic1};

%% Define the terms of the PDE
%
% Here we have 2 terms, having only nDims=2 (x,y) operators.

dV = @(x,p,t,dat) 0*x+1;

%%
% -df/dx which is 
%
% d/dx g1(x) f(x,y)          [grad,g1(x)=-1, BCL=P, BCR=P]

g1 = @(x,p,t,dat) x*0-1;
pterm1  =  DIV(num_dims,g1,'',-1,'P','P','','','',dV);
term1_x = SD_TERM({pterm1});

g1 = @(x,p,t,dat) 0*x+1;
pterm1 = MASS(g1,'','',dV);
term1_y = SD_TERM({pterm1});

term1   = MD_TERM(num_dims,{term1_x,term1_y});

%%
% -df/dy which is
%
% d/dy g1(y) f(x,y)          [grad,g1(y)=-1, BCL=P, BCR=P]

g1 = @(y,p,t,dat) y*0-1;
pterm1  =  DIV(num_dims,g1,'',-1,'P','P','','','',dV);
term2_y = SD_TERM({pterm1});
term2   = MD_TERM(num_dims,{[],term2_y});

terms = {term1,term2};

%% Define some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

%% Define sources

%%
% Source term 1
source1 = { ...
    @(x,p,t) cos(pi*x),   ...   % s1x
    @(y,p,t) sin(2*pi*y), ...   % s1y
    @(t)   2*cos(2*t)   ...   % s1t
    };

%%
% Source term 2
source2 = { ...
    @(x,p,t)  cos(pi*x),    ...   % s2x
    @(y,p,t)  cos(2*pi*y),  ...   % s2y
    @(t)    2*pi*sin(2*t) ...   % s2t
    };

%%
% Source term 3
source3 = { ...
    @(x,p,t)  sin(pi*x),   ...  % s3x
    @(y,p,t)  sin(2*pi*y), ...  % s3y
    @(t)    -pi*sin(2*t) ...  % s3t
    };

sources = {source1,source2,source3};

%% Define the analytic solution (optional).

soln1 = new_md_func(num_dims,{soln_x,soln_y,soln_t});
solutions = {soln1};

%% Define function to set dt

    function dt=set_dt(pde,CFL)
        
        Lmax = pde.dimensions{1}.max;
        Lmin = pde.dimensions{1}.min;
        LevX = pde.dimensions{1}.lev;
        dt = (Lmax-Lmin)/2^LevX*CFL;
    end

%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,[],initial_conditions,solutions);

end


