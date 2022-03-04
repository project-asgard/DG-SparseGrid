function pde = lb_collision2(opts)
% Lenard-Bernstein Equation in 2D
%
% PDE:
%
% df/dt = div( (v-u)f + th\grad f)
%
% Domain is [vx_min,vx_max]x[vy_min,vy_max]
% Dirichlet boundary condition in advection
% Neumann boundary condition in diffusion
%
% Diffusion terms are dealt with via LDG
%
% Run with
%
% explicit
% asgard(@lb_collision2D,'CFL',0.01);
%
% implicit
% asgard(@lb_collision2,'timestep_method','BE','dt',0.2,'num_steps',20,'case',2,'lev',4,'deg',3,'grid_type','SG')
%
% cyclotron motion test case
% asgard(@lb_collision2,'timestep_method','BE','dt',0.2,'num_steps',10,'case',4,'lev',3,'deg',6,'grid_type','SG','time_independent_build_A',true)
%
% case = 1 % Initial Condition: Square
% case = 2 % Initial Condition: Maxwellian
% case = 3 % Initial Condition: Double Maxwellian
% case = 4 % Cyclotron motion only

%% Macro-Micro?

MM = false;

%% Define the parameters

switch opts.case_
    case 1
        params.n  = 1.0; % Density
        params.ux = 0.0; % Velocity (x)
        params.uy = 0.0; % Velocity (y)
        params.th = 1.0/3.0; % Temperature
        vx_min = -4.0;
        vx_max = +4.0;
        vy_min = -4.0;
        vy_max = +4.0;
    case 2
        params.n  = 1.0; % Density
        params.ux = 1.0; % Velocity (x)
        params.uy = 1.0; % Velocity (y)
        params.th = 1.5; % Temperature
        vx_min = -5.0;
        vx_max = +5.0;
        vy_min = -5.0;
        vy_max = +5.0;
    case 3
        params.n  = 2.0; % Density
        params.ux = 1.0; % Velocity (x)
        params.uy = 1.0; % Velocity (y)
        params.th = 0.5; % Temperature
        vx_min = -5.0;
        vx_max = +5.0;
        vy_min = -5.0;
        vy_max = +5.0;
    case 4
        params.Bz = 2;
        params.n  = 1.0; % Density
        params.ux = 0.0; % Velocity (x)
        params.uy = 0.0; % Velocity (y)
        params.th = 1.0/3.0; % Temperature
        vx_min = -4.0;
        vx_max = +4.0;
        vy_min = -4.0;
        vy_max = +4.0;
    otherwise
end

%% Define the dimensions

dim_x = DIMENSION(vx_min,vx_max);
dim_x.moment_dV = @(x,p,t,dat) 0*x+1;
dim_y = DIMENSION(vy_min,vy_max);
dim_y.moment_dV = @(x,p,t,dat) 0*x+1;
dimensions = {dim_x,dim_y};
num_dims = numel(dimensions);

%% Define the analytic solution (optional).

switch opts.case_
    case 1
        if MM
            soln_x = @(x,p,t) 0*x;
            soln_y = @(y,p,t) 0*y;
        else
            soln_x = @(x,p,t) sqrt(p.n/(2*pi*p.th))*exp(-(x-p.ux).^2/(2*p.th));
            soln_y = @(y,p,t) sqrt(p.n/(2*pi*p.th))*exp(-(y-p.uy).^2/(2*p.th));
        end
        soln_t = @(t,p) 0*t+1;
        soln1  = new_md_func(num_dims,{soln_x,soln_y,soln_t});
    case 2
        soln_x = @(x,p,t) sqrt(p.n/(2*pi*p.th))*exp(-(x-p.ux).^2/(2*p.th));
        soln_y = @(y,p,t) sqrt(p.n/(2*pi*p.th))*exp(-(y-p.uy).^2/(2*p.th));
        soln_t = @(t,p) 0*t+1;
        soln1  = new_md_func(num_dims,{soln_x,soln_y,soln_t});
    case 3
        soln_x = @(x,p,t) sqrt(p.n/(2*pi*p.th))*exp(-(x-p.ux).^2/(2*p.th));
        soln_y = @(y,p,t) sqrt(p.n/(2*pi*p.th))*exp(-(y-p.uy).^2/(2*p.th));
        soln_t = @(t,p) 0*t+1;
        soln1  = new_md_func(num_dims,{soln_x,soln_y,soln_t});
    case 4
        if MM
            soln_x = @(x,p,t) 0*x;
            soln_y = @(y,p,t) 0*y;
        else
            soln_x = @(x,p,t) sqrt(p.n/(2*pi*p.th))*exp(-(x-p.ux).^2/(2*p.th));
            soln_y = @(y,p,t) sqrt(p.n/(2*pi*p.th))*exp(-(y-p.uy).^2/(2*p.th));
        end
        soln_t = @(t,p) 0*t+1;
        soln1  = new_md_func(num_dims,{soln_x,soln_y,soln_t});
    otherwise
end

solutions = {soln1};

%% Define the boundary conditions

BCL = new_md_func(num_dims,{@(x,p,t,dat) 0*x,@(y,p,t) 0*y,@(t,p) 0*t+1});
BCR = new_md_func(num_dims,{@(x,p,t,dat) 0*x,@(y,p,t) 0*y,@(t,p) 0*t+1});

%% Initial conditions

switch opts.case_
    case 1
        ic_x = @(x,p,t,d) 0.5.*((x > -1.0) - (x > +1.0));
        ic_y = @(y,p,t,d) 0.5.*((y > -1.0) - (y > +1.0));
        ic_t = @(t,p) 0*t+1;
        ic1 = new_md_func(num_dims,{ic_x,ic_y,ic_t});
        if MM
            ic_x = @(x,p,t,d) - sqrt(p.n/(2*pi*p.th))*exp(-(x-p.ux).^2/(2*p.th));
            ic_y = @(y,p,t,d)   sqrt(p.n/(2*pi*p.th))*exp(-(y-p.uy).^2/(2*p.th));
        else
            ic_x = @(x,p,t,d) 0*x;
            ic_y = @(y,p,t,d) 0*y;
        end
        ic2 = new_md_func(num_dims,{ic_x,ic_y,ic_t});
        initial_conditions = {ic1,ic2};
    case 2
        ic_x = soln_x;
        ic_y = soln_y;
        ic_t = soln_t;
        ic = new_md_func(num_dims,{ic_x,ic_y,ic_t});
        initial_conditions = {ic};
    case 3
        ic1_x = @(x,p,t) 1.0/sqrt(pi)*exp(-(x-1.0).^2);
        ic1_y = @(y,p,t) 1.0/sqrt(pi)*exp(-(y-0.0).^2);
        ic2_x = @(x,p,t) 1.0/sqrt(pi)*exp(-(x-0.0).^2);
        ic2_y = @(y,p,t) 1.0/sqrt(pi)*exp(-(y-1.0).^2);
        ic_t  = @(t,p) 0*t+1;
        ic1   = new_md_func(num_dims,{ic1_x,ic1_y,ic_t});
        ic2   = new_md_func(num_dims,{ic2_x,ic2_y,ic_t});
        initial_conditions = {ic1,ic2};
    case 4
        sig = 0.5;
        ic1_x = @(x,p,t) 1.0/sqrt(pi)*exp(-(x-1.0).^2/sig);
        ic1_y = @(y,p,t) 1.0/sqrt(pi)*exp(-(y-0.0).^2/sig);
        ic_t  = @(t,p) 0*t+1;
        ic1   = new_md_func(num_dims,{ic1_x,ic1_y,ic_t});
        initial_conditions = {ic1};
    otherwise
        ic_x = soln_x;
        ic_y = soln_y;
        ic_t = soln_t;
        ic = new_md_func(num_dims,{ic_x,ic_y,ic_t});
        initial_conditions = {ic};
end

%% Construct moments

% Mass Moment
moment_x = @(x,p,t) 0*x+1;
moment_y = @(y,p,t) 0*y+1;
moment_t = @(p,t)   0*t+1;
moment_func = new_md_func(num_dims,{moment_x,moment_y,moment_t});
moment0 = MOMENT({moment_func});

% Momentum Moment (x)

moment_x = @(x,p,t) x;
moment_y = @(y,p,t) 0*y+1;
moment_t = @(p,t)   0*t+1;
moment_func = new_md_func(num_dims,{moment_x,moment_y,moment_t});
moment1 = MOMENT({moment_func});

% Momentum Moment (y)

moment_x = @(x,p,t) 0*x+1;
moment_y = @(y,p,t) y;
moment_t = @(p,t)   0*t+1;
moment_func = new_md_func(num_dims,{moment_x,moment_y,moment_t});
moment2 = MOMENT({moment_func});

% Energy Moment

moment_x1 = @(x,p,t) x.^2;
moment_y1 = @(y,p,t) 0*y+1;
moment_x2 = @(x,p,t) 0*x+1;
moment_y2 = @(y,p,t) y.^2;
moment_t  = @(p,t)   0*t+1;
moment_func1 = new_md_func(num_dims,{moment_x1,moment_y1,moment_t});
moment_func2 = new_md_func(num_dims,{moment_x2,moment_y2,moment_t});
moment3 = MOMENT({moment_func1,moment_func2});

moments = {moment0,moment1,moment2,moment3};

%% Define the terms of the PDE
%

dV = @(x,p,t,dat) 0*x+1;

%%
% Setup the advection terms div( (v-u)f ) = d/dx((vx-ux)f) + d/dy((vy-uy)f)

%%
% d/dx((vx-ux)f)

g1 = @(x,p,t,dat) x-p.ux;
pterm1  = DIV(num_dims,g1,'',-1,'D','D',BCL,BCR,'',dV);
term1_x = SD_TERM({pterm1});
term1   = MD_TERM(num_dims,{term1_x,[]});

%%
% d/dy((vy-uy)f)
g1 = @(y,p,t,dat) y-p.uy;
pterm1  = DIV(num_dims,g1,'',-1,'D','D',BCL,BCR,'',dV);
term2_y = SD_TERM({pterm1});
term2   = MD_TERM(num_dims,{[],term2_y});

%%
% Setup the diffusion term div(theta\grad f) term
%
% div(theta\grad f) is the system
%
% div(sqrt(theta)q) = d/dx(sqrt(theta)qx)+d/dy(sqrt(theta)qy)
%
% qx = sqrt(theta)df/dx
%
% qy = sqrt(theta)df/dy
%
%%
% Set up df/dt = d/dx(sqrt(theta)qx)
%           qx = d/dx(sqrt(theta)f)

g1 = @(x,p,t,dat) x.*0+sqrt(p.th);
g2 = @(x,p,t,dat) x.*0+sqrt(p.th);

pterm1  =  DIV(num_dims,g1,'',+1,'D','D',BCL,BCR,'',dV);
pterm2  = GRAD(num_dims,g2,'',-1,'N','N','' ,'' ,'',dV);
term3_x = SD_TERM({pterm1,pterm2});

g1 = @(x,p,t,dat) x.*0+1.0;
pterm1  = MASS(g1);
term3_y = SD_TERM({pterm1,pterm1});

term3   = MD_TERM(num_dims,{term3_x,term3_y});

%%
% Set up df/dt = d/dy(sqrt(theta)qy)
%           qy = d/dy(sqrt(theta)f)

g1 = @(y,p,t,dat) y.*0+1.0;

pterm1  = MASS(g1);
term4_x = SD_TERM({pterm1,pterm1});

g1 = @(y,p,t,dat) y.*0+sqrt(p.th);
g2 = @(y,p,t,dat) y.*0+sqrt(p.th);

pterm1  =  DIV(num_dims,g1,'',+1,'D','D',BCL,BCR,'',dV);
pterm2  = GRAD(num_dims,g2,'',-1,'N','N','' ,'' ,'',dV);
term4_y = SD_TERM({pterm1,pterm2});

term4 = MD_TERM(num_dims,{term4_x,term4_y});

%%
% v x B transport (cyclotron motion around a prescribed B vector B=Bz*zhat
% DLG Note: we have to split each single dimension advection term up into 
%           2 terms to allow for upwinding of different signs, resulting in 4 terms
%           (term5a,term5b and term6a,term6b)

% -d/dvx(vy*Bz*f) => vy*Bz * d/dvx(f)

g = @(vx,p,t,d) vx.*0+1;
pterm = DIV(num_dims,g,'',-1,'N','D',BCL,BCR,'',dV);
term5_vx = SD_TERM({pterm});
g = @(vy,p,t,d) (-vy.*p.Bz) .* (vy<=0);
pterm = MASS(g);
term5_vy = SD_TERM({pterm});
term5a = MD_TERM(num_dims,{term5_vx,term5_vy});

g = @(vx,p,t,d) vx.*0+1;
pterm = DIV(num_dims,g,'',+1,'D','N',BCL,BCR,'',dV);
term5_vx = SD_TERM({pterm});
g = @(vy,p,t,d) (-vy.*p.Bz) .* (vy>0);
pterm = MASS(g);
term5_vy = SD_TERM({pterm});
term5b = MD_TERM(num_dims,{term5_vx,term5_vy});

% +d/dvy(vx*Bz*f) => vx*Bz * d/dvy(f)

g = @(vx,p,t,d) (vx.*p.Bz) .* (vx>=0);
pterm = MASS(g);
term6_vx = SD_TERM({pterm});
g = @(vy,p,t,d) vy.*0 + 1;
pterm = DIV(num_dims,g,'',-1,'N','D',BCL,BCR,'',dV);
term6_vy = SD_TERM({pterm});
term6a = MD_TERM(num_dims,{term6_vx,term6_vy});

g = @(vx,p,t,d) (vx.*p.Bz) .* (vx<0);
pterm = MASS(g);
term6_vx = SD_TERM({pterm});
g = @(vy,p,t,d) vy.*0 + 1;
pterm = DIV(num_dims,g,'',+1,'D','N',BCL,BCR,'',dV);
term6_vy = SD_TERM({pterm});
term6b = MD_TERM(num_dims,{term6_vx,term6_vy});

switch opts.case_
    case 4
        terms = {term5a,term5b,term6a,term6b};
    otherwise
        terms = {term1,term2,term3,term4};
end


%% Define some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

%% Define sources

sources = {};

%% Define a function to set dt

    function dt=set_dt(pde,CFL)
        dims = pde.dimensions;
        % for Diffusion equation: dt = C * dx^2
        lev = dims{1}.lev;
        dx = 1/2^lev;
        dt = CFL*dx^2;
    end

%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,[],initial_conditions,solutions,moments);

end