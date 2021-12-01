function pde = fokkerplanck2_C_div(opts)
% Combining momentum and pitch angle dynamics
%
% Problems 6.1, 6.2, and 6.3 from the RE paper.
%
% d/dt f(p,z,t) == div( A^2*F_C1 + F_C2 )
%
% where
%
% \F_C1 = grad(f) and \F_C2 = C_F(p)f\hat{p}
%
% Here 
%
%     div(F_p\hat{p} + F_z\hat{z}) = 1/p^2*d/dp(p^2F_p)
%                                               + 1/p*d/dz(sqrt(1-z^2)F_z)
% and
%
%     grad(f) = df/dp*hat{p} + 1/p*sqrt(1-z^2)*df/dz*hat{z}.
%
% The diffusion tensor A^2 is given by
%
% A\hat{p} = sqrt(C_A(p))\hat{p} and A\hat{z} = 1/p*sqrt(C_B(p))\hat{z}.
%
% MIXED formulation:
%
% d/dt f(p,z,t) == -div( A*q + Gamma_C2 )
%
% q = A*grad(f) = q_p\hat{p} + q_z\hat{z}
%
% Run with
%
% explicit
% asgard(@fokkerplanck2_C_div,'CFL',0.01,'case',1)
%
% implicit
% asgard(@fokkerplanck2_C_div,'timestep_method','BE','num_steps',20,'CFL',1.0,'deg',3,'lev',4,'case',1)
%
% case = 1 % step function in momentum
% case = 2 %
% case = 3 % 

params = fokkerplanck_parameters(opts);

%% Setup the dimensions

dim_p = DIMENSION(0,+10);
dim_p.moment_dV = @(x,p,t) x.^2;
dim_z = DIMENSION(-1,+1);
dim_z.moment_dV = @(x,p,t) 0*x+1;
dimensions = {dim_p,dim_z};
num_dims = numel(dimensions);

%% Define the analytic solution (optional)

soln_p = @(x,p,t) p.soln_p(x,p,t);
soln_z = @(x,p,t) p.soln_z(x,p,t);
soln1 = new_md_func(num_dims,{soln_p,soln_z});
solutions = {soln1};

%% Define the initial conditions

ic_p = @(x,p,t) p.f0_p(x);
ic_z = @(x,p,t) p.f0_z(x);
ic1 = new_md_func(num_dims,{ic_p,ic_z});
initial_conditions = {ic1};

%% Define the boundary conditions

%% Look at some of the coefficients for sanity checking

p = 0:.001:10;
plot(p,p.^2.*params.Ca(p));
hold on
plot(p,p.^2.*params.Cf(p));
plot(p,min(params.Cb(p)./p.^4,10));
hold off

%% Define the terms of the PDE

%%
% termC1 = LDG for q_p variable.
%
% LDG in q_p variable gives LDG in p, MASS in z
%
% becomes
%
% LDG in p
%       div( A*q_p*\hat{p} )     [div , g1(p) = sqrt(C_A), BCL=D,BRC=N]
%   q(p) == A*grad(f)*\hat{p}    [grad, g2(p) = sqrt(C_A), BCL=N,BCR=D]
%
% MASS in z (first eqn)             [mass, g1(z) = 1, BC=NA] 
%           (secnd eqn)             [mass, g1(z) = 1, BC=NA] 

% Surface jacobian is for p constant is p^2
dV_p = @(x,p,t,dat) x.^2;
dV_z = @(x,p,t,dat) 0*x+1;

% LDG in p
g1 = @(x,p,t,dat) sqrt(p.Ca(x));
pterm1  =  DIV(num_dims,g1,'',+1,'D','N','','','',dV_p);
pterm2  = GRAD(num_dims,g1,'',-1,'N','D','','','',dV_p);
term1_p = SD_TERM({pterm1,pterm2});

% MASS in z
g1 = @(x,p,t,dat) 0*x+1;
pterm1 = MASS(g1,'','',dV_z);
term1_z = SD_TERM({pterm1,pterm1});

termC1  = MD_TERM(num_dims,{term1_p,term1_z});

%%
% termC2 == -div(Gamma_C2) 
%
% Gamma_C2 = -C_F(p)f\hat{p}  <-- Flow is right to left
%
% This term is DIV in p and MASS in z
%
% DIV in p
%     == div( C_F(p)f\hat{p})      [div, g1(p)=C_f(p), BCL=N, BCR=D]
%
% MASS in z                        [mass, g1(z) = 1, BC=NA] 

% Surface jacobian for p constant is p^2
dV_p = @(x,p,t,dat) x.^2;
dV_z = @(x,p,t,dat) 0*x+1;

% DIV in p
g1 = @(x,p,t,dat) p.Cf(x);
pterm1  =  DIV(num_dims,g1,'',-1,'N','D','','','',dV_p);
term2_p = SD_TERM({pterm1});

% MASS in z
g1 = @(x,p,t,dat) 0*x+1;
pterm1 = MASS(g1,'','',dV_z);
term2_z = SD_TERM({pterm1});

termC2   = MD_TERM(num_dims,{term2_p,term2_z});

%%
% termC3 = LDG for q_z variable.
%
% LDG in q_z variable gives MASS in p, LDG in z
%
% becomes
% MASS in p  (first eqn)               [mass, g1(p) = sqrt(C_b(p)), BC=NA] 
%            (secnd eqn)               [mass, g1(p) = sqrt(C_b(p)), BC=NA] 
%
% LDG in z
%       div( A*q_z*\hat{z} )           [div , g1(p) = 1, BCL=D,BRC=N]
%   q(p) == A*grad(f)*\hat{z}          [grad, g2(p) = 1, BCL=N,BCR=D]
%

% Surface jacobian for z constnat is p*sqrt(1-z^2)
dV_p = @(x,p,t,dat) x;
dV_z = @(x,p,t,dat) sqrt(1-x.^2);

% MASS in p

g1 = @(x,p,t,dat) sqrt(p.Cb(x));
pterm1  = MASS(g1,'','',dV_p);
term3_p = SD_TERM({pterm1,pterm1});

% LDG in z
g1 = @(x,p,t,dat) x.*0+1;
pterm1  =  DIV(num_dims,g1,'',+1,'D','D','','','',dV_z);
pterm2  = GRAD(num_dims,g1,'',-1,'N','N','','','',dV_z);
term3_z = SD_TERM({pterm1,pterm2});

termC3 = MD_TERM(num_dims,{term3_p,term3_z});

terms = {termC1,termC2,termC3};

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

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,[],initial_conditions,solutions);

end


