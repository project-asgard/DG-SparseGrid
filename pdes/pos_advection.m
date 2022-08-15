function pde = pos_advection(opts)
% 2D test case using continuity equation, i.e.,
%
% df/dt == -df/dx - df/dy
%       == -div(cf) where c=(c1,c2)=(1,1)
%
% Run with
%
% explicit
% asgard(@pos_advection,'lev',6,'deg',1,'grid_type','FG','timestep_method','IMEX','dt',0.01,'num_steps',50,'quiet',false,'build_realspace_output',false,'plot_freq',10,'adapt_threshold',1e-3,'calculate_mass',false,'max_lev',6,'fast_FG_matrix_assembly',true)
% asgard(@pos_advection,'lev',6,'deg',1,'grid_type','SG','timestep_method','IMEX','dt',0.005,'num_steps',50,'quiet',false,'build_realspace_output',false,'plot_freq',10,'calculate_mass',false,'max_lev',6,'fast_FG_matrix_assembly',true)

% pde.CFL = 0.1;

%% Define the dimensions
%
% Here we setup a 2D problem (x,y)

dim_x = DIMENSION(-1,+1);
dim_x.moment_dV = @(x,p,t,dat) 0*x+1;
dim_y = DIMENSION(-1,+1);
dim_y.moment_dV = @(y,p,t,dat) 0*y+1;
dimensions = {dim_x,dim_y};
num_dims = numel(dimensions);

%% Solution

%soln_x = @(x,p,t) (x-t > 0).*(x-t < 0.25);
%soln_y = @(y,p,t) (y-t > -0.25).*(y-t < 0);
%  Box Support
soln_x = @(x,p,t) (keep_bound(x-t) > 0).*(keep_bound(x-t) < 0.5);
soln_y = @(y,p,t) (keep_bound(y-t) > -0.5).*(keep_bound(y-t) < 0);
%  
theta = 0.5;
%soln_x = @(x,p,t) exp(-(keep_bound(x-t)/theta).^2);
%soln_y = @(x,p,t) exp(-(keep_bound(x-t)/theta).^2);
%soln_x = @(x,p,t) 1-keep_bound(x-t).^2;
%soln_y = @(x,p,t) 1-keep_bound(x-t).^2;
%soln_x = @(x,p,t) x;
%soln_y = @(y,p,t) 0*y+1;
soln_t = @(t,p) 0*t+1;

%% Construct moments

%mass moment
moment_func = new_md_func(num_dims,{@(x,p,t) 0*x+1,...
                                    @(x,p,t) 0*x+1,...
                                    @(p,t)   0*t+1});
moment0 = MOMENT({moment_func});

moments = {moment0};

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

term1   = MD_TERM(num_dims,{term1_x,term1_y},'E');

%%
% -df/dy which is
%
% d/dy g1(y) f(x,y)          [grad,g1(y)=-1, BCL=P, BCR=P]

g1 = @(y,p,t,dat) y*0-1;
pterm1  =  DIV(num_dims,g1,'',-1,'P','P','','','',dV);
term2_y = SD_TERM({pterm1});

term2   = MD_TERM(num_dims,{[],term2_y},'E');

terms = {term1,term2};
%terms = {};

%% Define some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

%% Define sources

sources = {};

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

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,[],initial_conditions,solutions,moments);

end

function z = keep_bound(y)
    zeta = floor(y/2+0.5);
    z = y - 2*zeta;
end



