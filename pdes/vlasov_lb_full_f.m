function pde = vlasov_lb_full_f(opts)
% 2D test case using continuitv equation, i.e.,
%
% df/dt == -v*\grad_x f + div_v( (v-u)f + theta\grad_v f)
%
% BC is peridoic in x
% BC in v is all inflow in advection for v and Neumann for diffusion in v
%
% Run with
%
% explicit
% asgard(@vlasov_lb_full_f,'lev',3,'deg',3,'CFL',0.1)
%
% implicit
% asgard(@vlasov_lb_full_f,'timestep_method','BE','deg',3,'lev',3,'dt',0.1)
%
% with adaptivitv
% asgard(@vlasov_lb_full_f,'timestep_method','CN','adapt',true)

soln_x = @(x,p,t)  0*x+1;
soln_v = @(v,p,t)  0*v+1;
soln_t = @(t,p)    0*t+1;

%% Define the dimensions
%
% Here we setup a 2D problem (x,v)

dim_x = DIMENSION(-1,+1);
dim_x.moment_dV = @(x,p,t,dat) 0*x+1;
dim_v = DIMENSION(-6,+6);
dim_v.moment_dV = @(v,p,t,dat) 0*v+1;
dimensions = {dim_x,dim_v};
num_dims = numel(dimensions);

%% Construct moments

%mass moment
moment_func = new_md_func(num_dims,{@(x,p,t) 0*x+1,...
                                    @(x,p,t) 0*x+1,...
                                    @(p,t)   0*t+1});
moment0 = MOMENT({moment_func});

%momentum moment
moment_func = new_md_func(num_dims,{@(x,p,t) 0*x+1,...
                                    @(x,p,t) x,...
                                    @(p,t)   0*t+1});
moment1 = MOMENT({moment_func});

%energy moment
moment_func = new_md_func(num_dims,{@(x,p,t) 0*x+1,...
                                    @(x,p,t) x.^2,...
                                    @(p,t)   0*t+1});
moment2 = MOMENT({moment_func});

moments = {moment0,moment1,moment2};

%% Construct (n,u,theta)

params.n  = @(x) 1*(x<-0.5) + 1/8*(x >= -0.5).*(x <= 0.5) + 1*(x>0.5);
params.u  = @(x) 0;
params.th = @(x) 1*(x<-0.5) + 5/4*(x >= -0.5).*(x <= 0.5) + 1*(x>0.5);
params.nu = 1e-3;

%% Define the analvtic solution (optional).

%soln1 = new_md_func(num_dims,{soln_x,soln_v,soln_t});
%solutions = {soln1};

solutions = {};

%% Initial conditions
ic1 = new_md_func(num_dims,{@(x,p,t) (abs(x) > 0.5),...
                            @(v,p,t) 1/sqrt(2*pi)*exp(-v.^2/2),...
                            @(t,p) 0*t+1});
ic2 = new_md_func(num_dims,{@(x,p,t) (abs(x) <= 0.5),...
                            @(v,p,t) (1/8)/sqrt(2*pi*5/4)*exp(-v.^2/(2*5/4)),...
                            @(t,p) 0*t+1});                        
initial_conditions = {ic1,ic2};

%% Define the terms of the PDE
%
% Here we have 2 terms, having onlv nDims=2 (x,v) operators.

dV = @(x,p,t,dat) 0*x+1;

%% Term 1
% -v\cdot\grad_x f for v > 0
%

g1 = @(x,p,t,dat) x*0-1;
pterm1  =  DIV(num_dims,g1,'',-1,'P','P','','','',dV);
term_x = SD_TERM({pterm1});

g1 = @(x,p,t,dat) x.*(x>0);
pterm1 = MASS(g1,'','',dV);
term_v = SD_TERM({pterm1});

term1   = MD_TERM(num_dims,{term_x,term_v});

%% Term 2
% -v\cdot\grad_x f for v < 0
%

g1 = @(x,p,t,dat) x*0-1;
pterm1  =  DIV(num_dims,g1,'',+1,'P','P','','','',dV);
term_x = SD_TERM({pterm1});

g1 = @(x,p,t,dat) x.*(x<0);
pterm1 = MASS(g1,'','',dV);
term_v = SD_TERM({pterm1});

term2   = MD_TERM(num_dims,{term_x,term_v});

%% Term 3
% v\cdot\grad_v f
%

g1 = @(x,p,t,dat) x*0+1/p.nu;
pterm1 = MASS(g1,'','',dV);
term_x = SD_TERM({pterm1});

g1 = @(x,p,t,dat) x;
pterm1  =  DIV(num_dims,g1,'',-1,'D','D','','','',dV);
term_v = SD_TERM({pterm1});

term3   = MD_TERM(num_dims,{term_x,term_v});

%% Term 4
% -u\cdot\grad_v f
%

g1 = @(x,p,t,dat) -p.u(x);
pterm1 = MASS(g1,'','',dV);
term_x = SD_TERM({pterm1});

g1 = @(x,p,t,dat) 0*x+1/p.nu;
pterm1  =  DIV(num_dims,g1,'',0,'D','D','','','',dV);
term_v = SD_TERM({pterm1});

term4   = MD_TERM(num_dims,{term_x,term_v});

%% Term 5
% div_v(th\grad_v f)
%
% Split by LDG
%
% div_v(th q)
% q = \grad_v f

g1 = @(x,p,t,dat) 1/p.nu;
g2 = @(x,p,t,dat) 0*x+p.th(x);
pterm1 = MASS(g1,'','',dV);
pterm2 = MASS(g2,'','',dV);
term_x = SD_TERM({pterm1,pterm2});

g1 = @(x,p,t,dat) 0*x+1;
g2 = @(x,p,t,dat) 0*x+1;
pterm1  =  DIV(num_dims,g1,'',-1,'D','D','','','',dV);
pterm2  = GRAD(num_dims,g2,'',+1,'N','N','','','',dV);
term_v = SD_TERM({pterm1,pterm2});

term5   = MD_TERM(num_dims,{term_x,term_v});

%% Combine terms
terms = {term1,term2,term3,term4,term5};

%% Define sources

sources = {};

%% Define function to set dt

    function dt=set_dt(pde,CFL)
        
        Lmax = pde.dimensions{2}.max;
        Lmin = pde.dimensions{2}.min;
        LevX = pde.dimensions{2}.lev;
        dt = (Lmax-Lmin)/2^LevX*CFL;
    end

%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,[],initial_conditions,solutions,moments);

end


