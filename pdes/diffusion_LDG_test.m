function pde = diffusion_LDG_test(opts)

% Test 2D LDG diffusion in spherical coordinates
%
% df/dt == div(grad f) in spherical cordinates
%
% (p,th) is a spherical coordinate (p,th,ph), so the div and grad look like 
%
% div[A] = 1/p^2 * d/dp * p^2[A_p] + 1/(p*sin(th)) d/dth(sin(th) A_th
%
%      where A = A_p \hat{p} + A_th \hat{th}
%
% grad[f] = df/dp \hat{p} + 1/r df/dth \hat{th}
%
% and the volument_element dV = p^2*sin(th)
%
% Run with
%
% asgard(@diffusion_LDG_test,'timestep_method','BE','lev',3,'deg',4,'num_steps',50,'CFL',1.5)
%
% Analytic solution: f(r,th,t) = exp(-t)*(sin(r)/r^2-cos(r)/r)*cos(th)
%   using approprite Dirichlet BCs

%% Define some parameters and add to pde object.

params.parameter1 = 0;

%% Define the dimensions
dV_p = @(x,p,t,d) x.^2;
dim_p = DIMENSION(0,4.49340945790906); %Gives zero dirichlet BCs
%dim_p = DIMENSION(0.5,pi);
dim_p.moment_dV = dV_p;
dim_th = DIMENSION(0,pi);
%dim_th = DIMENSION(0.5,pi-0.5);
dim_th.moment_dV = @(x,p,t,d) sin(x);
dimensions = {dim_p,dim_th};
num_dims = numel(dimensions);

%% Define the analytic solution (optional).

%soln_p = @(x,p,t) soln(x);
soln1 = new_md_func(num_dims,{@(x,p,t) soln_p(x,p,t),...
                              @(x,p,t) soln_th(x,p,t),...
                              @(t,p) exp(-t) });
solutions = {soln1};

%% Define initial conditions

%ic_p = @(z,p,t) f0(z);
ic1 = new_md_func(num_dims,{@(x,p,t) soln_p(x,p,t),@(x,p,t) soln_th(x,p,t)});
initial_conditions = {ic1};

%% LHS terms (mass only)

g1 = @(x,p,t,dat) 0*x+1;
pterm1 = MASS(g1,[],[],dim_p.moment_dV);
pterm2 = MASS(g1,[],[],dim_th.moment_dV);
LHS_term_p = SD_TERM({pterm1});
LHS_term_th = SD_TERM({pterm2});

LHS_term = MD_TERM(num_dims,{LHS_term_p,LHS_term_th});

%LHS_terms = {LHS_term};
LHS_terms = {};

%% RHS terms

% term1 is done using SLDG in p with mass in th
%
% eq1 :  df/dt == div(q)        [pterm1: div (g1(p)=1,+1, BCL=N, BCR=N)]
% eq2 :      q == grad(f)       [pterm2: grad(g2(p)=1,-1, BCL=D, BCR=D)]
%

g1 = @(x,p,t,dat) 0*x+1;
dV_p = @(x,p,t,d) x.^2;
dV_th = @(x,p,t,d) sin(x);

pterm1 =  DIV(num_dims,g1,'',+1,'N','N','','','',dV_p);
pterm2 = GRAD(num_dims,g1,'',-1,'D','D',soln1,soln1,'',dV_p);
%pterm1 =  DIV(num_dims,g1,'',+1,'D','D','','','',dV_p);
%pterm2 = GRAD(num_dims,g1,'',-1,'N','N','','','',dV_p);
term1_p = SD_TERM({pterm1,pterm2}); % order here is as in equation

pterm1 = MASS(g1,[],[],dV_th); 
term1_th = SD_TERM({pterm1,pterm1});

term1   = MD_TERM(num_dims,{term1_p,term1_th});

% term2 is done using mass in p with SLDG in th.
%
% eq1 :  df/dt == div(q)        [pterm1: div (g1(p)=1,+1, BCL=N, BCR=N)]
% eq2 :      q == grad(f)       [pterm2: grad(g2(p)=1,-1, BCL=D, BCR=D)]
%

dV_p = @(x,p,t,d) x;
dV_th = @(x,p,t,d) sin(x);

pterm1 = MASS(g1,[],[],dV_p);
term2_p = SD_TERM({pterm1,pterm1});

pterm1 =  DIV(num_dims,g1,'',+1,'N','N','','','',dV_th);
pterm2 = GRAD(num_dims,g1,'',-1,'D','D',soln1,soln1,'',dV_th);
term2_th = SD_TERM({pterm1,pterm2});

term2 = MD_TERM(num_dims,{term2_p,term2_th});

terms = {term1,term2};

%% Define sources

sources = {};

%% Define function to set time step
    function dt=set_dt(pde,CFL)       
        dims = pde.dimensions;
        xRange = dims{1}.max-dims{1}.min;
        lev = dims{1}.lev;
        dx = xRange/2^lev;
        dt = CFL * dx;       
    end

%% Construct PDE

pde = PDE(opts,dimensions,terms,LHS_terms,sources,params,@set_dt,[],initial_conditions,solutions);

end

function z = soln_p(x,p,t)
    if abs(x) < 1e-7
        z = x/3 - x.^3/30 + x.^5/840; %Taylor expansion truncation near 0
    else
        z = sin(x)./(x.^2) - cos(x)./x;
    end
    %z = (x > pi/2-0.25).*(x < pi/2+0.25);
end

function z = soln_th(x,p,t)
    z = cos(x);
    %z = (x > pi/4-0.25).*(x < pi/4+0.25);
end


