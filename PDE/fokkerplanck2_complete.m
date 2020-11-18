function pde = fokkerplanck2_complete(opts)
% Combining momentum and pitch angle dynamics
%
% Full PDE from the 2D runaway electron paper 
%
% d/dt f(p,z) == -div(flux_C + flux_E + flux_R)
%
% where 
%
% flux_C is flux due to electron-ion collisions
% flux_E is the flux due to E accleration
% flux_R is the flux due to radiation damping
%
% -div(flux_C) == termC1 + termC2 + termC3
% 
% termC1 == 1/p^2*d/dp*p^2*Ca*df/dp
% termC2 == 1/p^2*d/dp*p^2*Cf*f
% termC3 == Cb(p)/p^4 * d/dz( (1-z^2) * df/dz )
%
% -div(flux_E) == termE1 + termE2
%
% termE1 == -E*z*f(z) * 1/p^2 (d/dp p^2 f(p))
% termE2 == -E*p*f(p) * d/dz (1-z^2) f(z)
%
% -div(flux_R) == termR1 + termR2
%
% termR1 == 1/p^2 d/dp p^2 gamma(p) p / tau f(p) * (1-z^2) * f(z)
% termR2 == -1/(tau*gam(p)) f(p) * d/dz z(1-z^2) f(z)
%
% Run with
%
% explicit
% asgard(@fokkerplanck2_complete,'CFL',0.01,'case',1)
%
% implicit
% asgard(@fokkerplanck2_complete,'timestep_method','CN','num_steps',20,'CFL',1.0,'deg',3,'lev',4,'case',1)
%
% with adaptivity
% asgard(@fokkerplanck2_complete,'timestep_method','CN','num_steps',20,'CFL',1.0,'deg',3,'lev',4, 'adapt', true,'case',1)
%
% NOTES
%
% For Lin's case the command that should reproduce her results
% asgard(@fokkerplanck2_complete,'timestep_method','BE','lev',6,'deg',5,'dt',0.01,'num_steps',3000,'time_independent_A',true,'case',1)
%
% David's adjustment can be run as follows ...
% asgard(@fokkerplanck2_complete,'timestep_method','BE','lev',3,'deg',6,'dt',1,'num_steps',50,'grid_type','FG','time_independent_A',true,'case',1)
% or, run with the high order time integrator ...
% asgard(f@okkerplanck2_complete,'timestep_method','ode15s','lev',3,'deg',5,'dt',50,'num_steps',1,'grid_type','SG','case',1)

params = fokkerplanck_parameters(opts);

%% Setup the dimensions 

dim_p = DIMENSION(0.1,+10);
dim_z = DIMENSION(-1,+1);
dim_p.jacobian = @(x,p,t) x.^2;
dim_z.jacobian = @(x,p,t) x.*0+1;
dimensions = {dim_p,dim_z};
num_dims = numel(dimensions);

%% Define the analytic solution (optional)

solutions = {};

%% Define the initial conditions
ic_p = @(x,p,t) p.f0_p(x);
ic_z = @(x,p,t) p.f0_z(x);
ic1 = new_md_func(num_dims,{ic_p,ic_z});
initial_conditions = {ic1};

%% Define the terms of the PDE

%% -div(flux_C) == termC1 + termC2 + termC3

%% 
% termC1 == 1/p^2*d/dp*p^2*Ca*df/dp
%
% becomes 
%
% termC1 == g1(p) q(p)        [mass, g1(p) = 1/p^2,  BC N/A]
%   q(p) == d/dp g2(p) r(p)   [grad, g2(p) = p^2*Ca, BCL=D,BCR=N]
%   r(p) == d/dp g3(p) f(p)   [grad, g3(p) = 1,      BCL=N,BCR=D]

g1 = @(x,p,t,dat) 1./x.^2;
g2 = @(x,p,t,dat) x.^2.*p.Ca(x);
g3 = @(x,p,t,dat) x.*0+1; 

pterm1  = MASS(g1);
pterm2  = GRAD(num_dims,g2,+1,'D','D');
pterm3  = GRAD(num_dims,g3,-1,'N','N');

term1_p = SD_TERM({pterm1,pterm2,pterm3});
termC1  = MD_TERM(num_dims,{term1_p,[]});

%%
% termC2 == 1/p^2*d/dp*p^2*Cf*f
%
% becomes
%
% termC2 == g1(p) q(p)       [mass, g1(p)=1/p^2,  BC N/A]
%   q(p) == d/dp g2(p) f(p)  [grad, g2(p)=p^2*Cf, BCL=N,BCR=D]

g1 = @(x,p,t,dat) 1./x.^2;
g2 = @(x,p,t,dat) x.^2.*p.Cf(x);

pterm1  = MASS(g1);
pterm2  = GRAD(num_dims,g2,-1,'N','N');

term2_p = SD_TERM({pterm1,pterm2});
termC2   = MD_TERM(num_dims,{term2_p,[]});

%%
% termC3 == Cb(p)/p^4 * d/dz( (1-z^2) * df/dz )
%
% becomes
%
% termC3 == q(p) r(z)
%   q(p) == g1(p)            [mass, g1(p) = Cb(p)/p^4, BC N/A]
%   r(z) == d/dz g2(z) s(z)  [grad, g2(z) = 1-z^2,     BCL=D,BCR=D]
%   s(z) == d/dz g3(z) f(z)  [grad, g3(z) = 1,         BCL=N,BCR=N]

g1 = @(x,p,t,dat) p.Cb(x)./x.^4;
pterm1  = MASS(g1);
term3_p = SD_TERM({pterm1});

g2 = @(x,p,t,dat) (1-x.^2);
g3 = @(x,p,t,dat) x.*0+1;
pterm1  = GRAD(num_dims,g2,+1,'D','D');
pterm2  = GRAD(num_dims,g3,-1,'N','N');

term3_z = SD_TERM({pterm1,pterm2});

termC3 = MD_TERM(num_dims,{term3_p,term3_z});

%% -div(flux_E) == termE1 + termE2

% termE1 == -E*z*f(z) * 1/p^2 (d/dp p^2 f(p))
%        == r(z) * q(p)
%   r(z) == g1(z) f(z)       [mass, g1(z) = -E*z,  BC N/A]
%   q(p) == g2(p) u(p)       [mass, g2(p) = 1/p^2, BC N/A]
%   u(p) == d/dp g3(p) f(p)  [grad, g3(p) = p^2,   BCL=N,BCR=D]

% g1 = @(x,p,t,dat) -E.*x;
% g2 = @(x,p,t,dat) 1./x.^2;
% g3 = @(x,p,t,dat) x.^2;

g1 = @(x,p,t,dat) -x;
g2 = @(x,p,t,dat) 1./x.^2;
g3 = @(x,p,t,dat) p.E*x.^2;

pterm1   = MASS(g1);
termE1_z = SD_TERM({pterm1});

pterm1 = MASS(g2); 
pterm2 = GRAD(num_dims,g3,0,'N','N');% Lin's Setting

termE1_p = SD_TERM({pterm1,pterm2});

termE1 = MD_TERM(num_dims,{termE1_p,termE1_z});

% termE2 == -E*p*f(p) * d/dz (1-z^2) f(z)
%        == q(p) * r(z)
%   q(p) == g1(p) f(p)       [mass, g1(p) = -E*p,  BC N/A]
%   r(z) == d/dz g2(z) f(z)  [grad, g2(z) = 1-z^2, BCL=N,BCR=N]

% g1 = @(x,p,t,dat) -E.*x;
% g2 = @(x,p,t,dat) 1-x.^2;
g1 = @(x,p,t,dat) -1./x;
g2 = @(x,p,t,dat) p.E*(1-x.^2);

pterm1   = MASS(g1);
termE2_p = SD_TERM({pterm1});

pterm1   = GRAD(num_dims,g2,+1,'N','N');% Lin's Setting

termE2_z = SD_TERM({pterm1});

termE2 = MD_TERM(num_dims,{termE2_p,termE2_z});

%% -div(flux_R) == termR1 + termR2

% termR1 == 1/p^2 d/dp p^2 gamma(p) p / tau f(p) * (1-z^2) * f(z)
%        == q(p) * r(z)
%   q(p) == g1(p) u(p)       [mass, g1(p) = 1/p^2,                BC N/A]
%   u(p) == d/dp g2(p) f(p)  [grad, g2(p) = p^3 * gamma(p) / tau, BCL=N,BCR=D]
%   r(z) == g3(z) f(z)       [mass, g3(z) = 1-z^2,                BC N/A]

g1 = @(x,p,t,dat) 1./x.^2;
g2 = @(x,p,t,dat) x.^3 .* p.gamma(x) ./ p.tau;
g3 = @(x,p,t,dat) 1-x.^2;

pterm1   = MASS(g1);% This is not needed - by Lin
pterm2   = GRAD(num_dims,g2,1,'N','N');% Lin's Setting

termR1_p = SD_TERM({pterm1,pterm2});

pterm1   = MASS(g3);
termR1_z = SD_TERM({pterm1});

termR1   = MD_TERM(num_dims,{termR1_p,termR1_z});

% termR2 == -1/(tau*gam(p)) f(p) * d/dz z(1-z^2) f(z)
%        == q(p) * r(z)
%   q(p) == g1(p) f(p)       [mass, g1(p) = -1/(tau*gamma(p)),    BC N/A]
%   r(z) == d/dz g2(z) f(z)  [grad, g2(z) = z(1-z^2),             BCL=N,BCR=N]

g1 = @(x,p,t,dat) -1./(p.tau.*p.gamma(x));
g2 = @(x,p,t,dat) x.*(1-x.^2);

pterm1   = MASS(g1);
termR2_p = SD_TERM({pterm1});

pterm1   = GRAD(num_dims,g2,0,'N','N');% Lin's Setting

termR2_z = SD_TERM({pterm1});

termR2 = MD_TERM(num_dims,{termR2_p, termR2_z});

terms = {termC1, termC2, termC3, termE1, termE2, termR1, termR2};


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


