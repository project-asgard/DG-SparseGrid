function pde = fokkerplanck2_C_LHS(opts)
% Combining momentum and pitch angle dynamics
%
% Problems 6.1, 6.2, and 6.3 from the RE paper.
%
% p^2 * d/dt f(p,z,t) == termC1 + termC2 + termC3
%
% termC1 == d/dp*p^2*Ca*df/dp
% termC2 == d/dp*p^2*Cf*f
% termC2 == Cb(p)/p^2 * d/dz( (1-z^2) * df/dz )
%
% Run with
%
% explicit
% asgard(@fokkerplanck2_C_LHS,'CFL',0.01,'case',1)
%
% implicit
% asgard(@fokkerplanck2_C_LHS,'timestep_method','CN','num_steps',20,'CFL',1.0,'deg',3,'lev',4,'case',1)
%
% NOTES
% DLG - case 3 doesn't seem to be working, although it does work for the
% nonLHS version of this PDE?

params = fokkerplanck_parameters(opts);

%% Define the dimensions 

dim_p = DIMENSION(0.01,+10);
dim_z = DIMENSION(-1,+1);
dimensions = {dim_p,dim_z};
num_dims = numel(dimensions);

%% Define the analytic solution (optional)

soln_p = @(x,p,t) p.soln_p(x,t);
soln_z = @(z,p,t) p.soln_z(z,t);
soln1 = new_md_func(num_dims,{soln_p,soln_z});
solutions = {soln1};

%% Define the initial conditions

ic_p = @(x,p,t) p.f0_p(x);
ic_z = @(x,p,t) p.f0_z(x);
ic1 = new_md_func(num_dims,{ic_p,ic_z});
initial_conditions = {ic1};

%% Define the terms of the PDE

%%
% LHS_term == p^2 d/dt f(p,z,t)

g1 = @(x,p,t,dat) x.^2;
pterm1 = MASS(g1);

LHS_term_p = SD_TERM({pterm1});
LHS_term   = MD_TERM(num_dims,{LHS_term_p,[]});

LHS_terms = {LHS_term};

%% 
% termC1 == d/dp*p^2*Ca*df/dp
%
% becomes 
%
% termC1 == d/dp g2(p) r(p)   [grad, g2(p) = p^2*Ca, BCL=D,BCR=N]        
%   r(p) == d/dp g3(p) f(p)   [grad, g3(p) = 1,      BCL=N,BCR=D]

g2 = @(x,p,t,dat) x.^2.*p.Ca(x);
g3 = @(x,p,t,dat) x.*0+1; 

pterm2  = GRAD(num_dims,g2,+1,'D','N');
pterm3  = GRAD(num_dims,g3,-1,'N','D');
term1_p = SD_TERM({pterm2,pterm3});
termC1  = MD_TERM(num_dims,{term1_p,[]});

%%
% termC2 == d/dp*p^2*Cf*f
%
% becomes
%
% termC2 == d/dp g2(p) f(p)  [grad, g2(p)=p^2*Cf, BCL=N,BCR=D]

g2 = @(x,p,t,dat) x.^2.*p.Cf(x);

pterm2  = GRAD(num_dims,g2,+1,'N','D');
term2_p = SD_TERM({pterm2});
termC2   = MD_TERM(num_dims,{term2_p,[]});

%%
% termC3 == Cb(p)/p^2 * d/dz( (1-z^2) * df/dz )
%
% becomes
%
% termC3 == q(p) r(z)
%   q(p) == g1(p)            [mass, g1(p) = Cb(p)/p^2, BC N/A]
%   r(z) == d/dz g2(z) s(z)  [grad, g2(z) = 1-z^2,     BCL=D,BCR=D]
%   s(z) == d/dz g3(z) f(z)  [grad, g3(z) = 1,         BCL=N,BCR=N]

g1 = @(x,p,t,dat) p.Cb(x)./x.^2;
pterm1  = MASS(g1);
term3_p = SD_TERM({pterm1});

g2 = @(x,p,t,dat) (1-x.^2);
g3 = @(x,p,t,dat) x.*0+1;
pterm1  = GRAD(num_dims,g2,+1,'D','D');
pterm2  = GRAD(num_dims,g3,-1,'N','N');
term3_z = SD_TERM({pterm1,pterm2});

termC3 = MD_TERM(num_dims,{term3_p,term3_z});

terms = {termC1,termC2,termC3};

%% Defin some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

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


