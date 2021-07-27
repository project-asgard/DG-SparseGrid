function pde = diffusion_LDG_test(opts)

% Test the momentum dynamics for the RE problem for E = R = 0
%
% df/dt == div(grad f) in spherical cordinates
%
% p is a spherical coordinate (p,th,ph), so the div and grad look like 
%
% div[] = 1/p^2 * d/dp * p^2[], grad[] = d/dp[]
%
% and the volument_element dV = p^2
%
% 
% split into two div(flux) terms (term1 and term2)
%
% term1 is done using SLDG defining eta1(p)=sqrt(psi(p)/p)
%
% eq1 :  df/dt == div(eta(p) * q)        [pterm1: div (g(p)=eta(p),+1, BCL=?, BCR=?)]
% eq2 :      q == eta(p) * grad(f)       [pterm2: grad(g(p)=eta(p),-1, BCL=D, BCR=N)]
%
% coeff_mat = pterm1.mat * pterm2.mat
%
% term2 is just a div
%
% eq1 :  df/dt == div(2*psi(p) * f)       [pterm1: div(g(p)=2*psi(p),+1, BCL=?, BCR=?]
%
% coeff_mat = pterm1.mat1
%
% Run with
%
% asgard(@diffusion_LDG_test,'timestep_method','BE','lev',3,'deg',4,'num_steps',50,'CFL',1.5,'case',1)
%
% case = 1 % maxwellian initial condition
% case = 2 % step function initial condition

%% Define some parameters and add to pde object.

params.parameter1 = 0;

%% Define the dimensions

dV_p = @(x,p,t,d) x.^2;
dim_p = DIMENSION(0,4.4934);
dim_p.moment_dV = dV_p;
dim_z = DIMENSION(0,pi/2);
dim_z.moment_dV = @(x,p,t,d) sin(x);
dimensions = {dim_p,dim_z};
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

LHS_terms = {};

%% RHS terms

% term1 is done using SLDG defining eta(p)=sqrt(psi(p)/p)
%
% eq1 :  df/dt == div(q)        [pterm1: div (g1(p)=1,+1, BCL=D, BCR=D)]
% eq2 :      q == grad(f)       [pterm2: grad(g2(p)=1,-1, BCL=N, BCR=N)]
%

g1 = @(x,p,t,dat) 0*x+1;

pterm1 =  DIV(num_dims,g1,'',+1,'N','N','','','',dV_p);
pterm2 = GRAD(num_dims,g1,'',-1,'D','D',soln1,soln1,'',dV_p);
term1_p = SD_TERM({pterm1,pterm2}); % order here is as in equation

pterm1 = MASS(g1,[],[],@(x,p,t,d) sin(x)); 
term1_th = SD_TERM({pterm1});
term1   = MD_TERM(num_dims,{term1_p,term1_th});

% term2 is just a div
%
% eq1 :  df/dt == div(2*psi(p) * f)       [pterm1: div(g3(p)=2*psi(p),+1, BCL=N, BCR=D]

dV_p = @(x,p,t,d) x;
dV_th = @(x,p,t,d) sin(x);
pterm1 = MASS(g1,[],[],dV_p);
term2_p = SD_TERM({pterm1});

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
end

function z = soln_th(x,p,t)
    z = cos(x);
end


