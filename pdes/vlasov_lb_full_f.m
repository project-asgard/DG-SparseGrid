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
%
% implicit
% asgard(@vlasov_lb_full_f,'timestep_method','IMEX','deg',3,'lev',[8 3],'dt',0.0002,'num_steps',500,'grid_type','FG','output_grid','interp','quiet',true,'build_realspace_output',false)
%

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

switch opts.case_ 
    case 1
        params.n  = @(x) 1*(x<-0.5) + 1/8*(x >= -0.5).*(x <= 0.5) + 1*(x>0.5);
        params.u  = @(x) 0;
        params.th = @(x) 1*(x<-0.5) + 4/5*(x >= -0.5).*(x <= 0.5) + 1*(x>0.5);
        params.nu = 1e3;
    case 2
        % deg = 3, lev = [8,5]
        U_0 = 1;
        T_0 = 1e-0;
        
        params.n  = @(x) 0.5*sin(2*pi*x)+1;
        params.u  = @(x) 0*x+U_0;
        params.th = @(x) 0*x+T_0;
        params.nu = 0;
    case 3
        N_0 = 1;
        T_0 = 1/3;
        C_0 = sqrt(3*T_0);
        Amp = 1e-6;
        
        params.n  = @(x) N_0 + Amp*sin(2*pi*x)/C_0^2;
        params.u  = @(x) Amp*sin(2*pi*x)./(C_0*params.n(x));
        params.th = @(x) (N_0*T_0+Amp*sin(2*pi*x))./params.n(x)-params.u(x).^2;
end

%% Define the analvtic solution (optional).

%soln1 = new_md_func(num_dims,{soln_x,soln_v,soln_t});
%solutions = {soln1};

switch opts.case_
    case 2 
        if params.nu < 1e-8
            soln1 = new_md_func(num_dims,...
                {@(x,p,t,dat) 0.5*sin(2*pi*x),...
                 @(v,p,t,dat) cos(2*pi*v*t).*1/sqrt(2*pi*T_0).*exp(-(v-U_0).^2/(2*T_0)),...
                 @(t,p) 0*t+1});
            soln2 = new_md_func(num_dims,...
                {@(x,p,t,dat) -0.5*cos(2*pi*x),...
                 @(v,p,t,dat) sin(2*pi*v*t).*1/sqrt(2*pi*T_0).*exp(-(v-U_0).^2/(2*T_0)),...
                 @(t,p) 0*t+1});
            soln3 = new_md_func(num_dims,...
                {@(x,p,t,dat) 0*x+1,...
                 @(v,p,t,dat) 1/sqrt(2*pi*T_0).*exp(-(v-U_0).^2/(2*T_0)),...
                 @(t,p) 0*t+1});
            solutions = {soln1,soln2,soln3};
        else
            solutions = {};
        end
    otherwise
        solutions = {};
end

%% Initial conditions

switch opts.case_
    case 1
        ic1 = new_md_func(num_dims,{@(x,p,t) (abs(x) > 0.5),...
                                   @(v,p,t) 1/sqrt(2*pi)*exp(-v.^2/2),...
                                   @(t,p) 0*t+1});
        ic2 = new_md_func(num_dims,{@(x,p,t) (abs(x) <= 0.5),...
                                   @(v,p,t) (1/8)/sqrt(2*pi*4/5)*exp(-v.^2/(2*4/5)),...
                                   @(t,p) 0*t+1});                        
        initial_conditions = {ic1,ic2};
    case 2 
        %U_0 and T_0 defined above
        ic1 = new_md_func(num_dims,{@(x,p,t) (0.5*sin(2*pi*x)+1)./sqrt(2*pi*T_0),...
                                   @(v,p,t) exp(-(v-U_0).^2/(2*T_0)),...
                                   @(t,p) 0*t+1});
        initial_conditions = {ic1};
end

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

term1   = MD_TERM(num_dims,{term_x,term_v},'E');

%% Term 2
% -v\cdot\grad_x f for v < 0
%

g1 = @(x,p,t,dat) x*0-1;
pterm1  =  DIV(num_dims,g1,'',+1,'P','P','','','',dV);
term_x = SD_TERM({pterm1});

g1 = @(x,p,t,dat) x.*(x<0);
pterm1 = MASS(g1,'','',dV);
term_v = SD_TERM({pterm1});

term2   = MD_TERM(num_dims,{term_x,term_v},'E');

%% Term 3
% v\cdot\grad_v f
%

g1 = @(x,p,t,dat) x*0+p.nu;
pterm1 = MASS(g1,'','',dV);
term_x = SD_TERM({pterm1});

g1 = @(x,p,t,dat) x;
pterm1  =  DIV(num_dims,g1,'',-1,'D','D','','','',dV);
term_v = SD_TERM({pterm1});

term3   = MD_TERM(num_dims,{term_x,term_v},'I');

%% Term 4
% -u\cdot\grad_v f
%

g1 = @(x,p,t,dat) -p.u(x);
pterm1 = MASS(g1,'','',dV);
term_x = SD_TERM({pterm1});

g1 = @(x,p,t,dat) 0*x+p.nu;
pterm1  =  DIV(num_dims,g1,'',0,'D','D','','','',dV);
term_v = SD_TERM({pterm1});

term4   = MD_TERM(num_dims,{term_x,term_v},'I');

%% Term 5
% div_v(th\grad_v f)
%
% Split by LDG
%
% div_v(th q)
% q = \grad_v f

g1 = @(x,p,t,dat) 0*x+1;
g2 = @(x,p,t,dat) 0*x+p.th(x)*p.nu;
pterm1 = MASS(g1,'','',dV);
pterm2 = MASS(g2,'','',dV);
term_x = SD_TERM({pterm1,pterm2});

g1 = @(x,p,t,dat) 0*x+1;
g2 = @(x,p,t,dat) 0*x+1;
pterm1  =  DIV(num_dims,g1,'',0,'D','D','','','',dV);
pterm2  = GRAD(num_dims,g2,'',0,'D','D','','','',dV);
term_v = SD_TERM({pterm1,pterm2});

term5   = MD_TERM(num_dims,{term_x,term_v},'I');

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


