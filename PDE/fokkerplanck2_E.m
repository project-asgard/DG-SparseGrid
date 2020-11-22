function pde = fokkerplanck2_E(opts)
% Combining momentum and pitch angle dynamics for the E term
%
% d/dt f(p,z) == -div(flux_E)
%
% flux_E is the flux due to E accleration
%
% -div(flux_E) == termE1 + termE2
%
% termE1 == -E*z*f(z) * 1/p^2 (d/dp p^2 f(p))
% termE2 == -E/p*f(p) * d/dz (1-z^2) f(z)
%
% Run with
%
% asgard(@fokkerplanck2_E,'timestep_method','matrix_exponential','dt',0.1,'num_steps',5,'case',2,'lev',4,'deg',4)

params = fokkerplanck_parameters(opts);

%% Setup the dimensions

dim_p = DIMENSION(0.1,+1);
dim_z = DIMENSION(-1,+1);
dim_p.jacobian = @(x,p,t) x.^2;
dim_z.jacobian = @(x,p,t) x.*0+1;
dimensions = {dim_p,dim_z};
num_dims = numel(dimensions);

%% Define the analytic solution (optional)

solutions = {};

%% Define the initial conditions

switch opts.case_
    case 1
        ic_p = @(x,p,t) x.*0+1;
        ic_z = @(z,p,t) z.*0+1;
    case 2
        ic_p = @(x,p,t) exp(-x.^2);
        ic_z = @(z,p,t) cos(z*pi)+1;
end

ic1 = new_md_func(num_dims,{ic_p,ic_z});
initial_conditions = {ic1};

%% Define the terms of the PDE

%% -div(flux_E) == termE1 + termE2

% termE1 == -E*z*f(z) * 1/p^2 (d/dp p^2 f(p))
%        == r(z) * q(p)
%   r(z) == g1(z) f(z)       [mass, g1(z) = -E*z,  BC N/A]
%   q(p) == g2(p) u(p)       [mass, g2(p) = 1/p^2, BC N/A]
%   u(p) == d/dp g3(p) f(p)  [grad, g3(p) = p^2,   BCL=N,BCR=D]

g1 = @(z,p,t,dat) -p.E .* z;
g2 = @(x,p,t,dat) min(1./x.^2,1e6);
g3 = @(x,p,t,dat) x.^2;
pterm1   = MASS(g1);
termE1_z = SD_TERM({pterm1});
pterm1   = MASS(g2);
pterm2   = GRAD(num_dims,g3,0,'N','N');% Lin's Setting (DLG-why is this central flux?)
termE1_p = SD_TERM({pterm1,pterm2});
termE1 = MD_TERM(num_dims,{termE1_p,termE1_z});

p = 0:.001:10;
plot(p,g1(p,params,[],[]));
hold on
plot(p,g2(p,params,[],[]));
plot(p,g3(p,params,[],[]));
hold off

% termE2 == -E/p*f(p) * d/dz (1-z^2) f(z)
%        == q(p) * r(z)
%   q(p) == g1(p) f(p)       [mass, g1(p) = -E*p,  BC N/A]
%   r(z) == d/dz g2(z) f(z)  [grad, g2(z) = 1-z^2, BCL=N,BCR=N]

g1 = @(x,p,t,dat) -p.E ./ x;
g2 = @(z,p,t,dat) 1-z.^2;
pterm1   = MASS(g1);
termE2_p = SD_TERM({pterm1});
pterm1   = GRAD(num_dims,g2,0,'N','N');% Lin's Setting
termE2_z = SD_TERM({pterm1});
termE2= MD_TERM(num_dims,{termE2_p,termE2_z});

% % termEA == -E/p df/dz
% g1 = @(x,p,t,dat) -p.E ./ x;
% pterm1 = MASS(g1);
% termEA_p = SD_TERM({pterm1});
% g1 = @(z,p,t,dat) z.*0+1;
% pterm1 = GRAD(num_dims,g1,+1,'N','N');
% termEA_z = SD_TERM({pterm1});
% termEA = MD_TERM(num_dims,{termEA_p,termEA_z});
% 
% % termEB == +E*z^2/p df/fz
% g1 = @(x,p,t,dat) +p.E ./ x;
% pterm1 = MASS(g1);
% termEB_p = SD_TERM({pterm1});
% g1 = @(z,p,t,dat) z.^2;
% g2 = @(z,p,t,dat) z.*0+1;
% pterm1 = MASS(g1);
% pterm2 = GRAD(num_dims,g2,+1,'N','N');
% termEB_z = SD_TERM({pterm1,pterm2});
% termEB = MD_TERM(num_dims,{termEB_p,termEB_z});
% 
% % termEC == +2E*z/p df/fz
% g1 = @(x,p,t,dat) +2*p.E ./ x;
% pterm1 = MASS(g1);
% termEC_p = SD_TERM({pterm1});
% g1 = @(z,p,t,dat) z;
% g2 = @(z,p,t,dat) z.*0+1;
% pterm1 = MASS(g1);
% pterm2 = GRAD(num_dims,g2,+1,'N','N');
% termEC_z = SD_TERM({pterm1,pterm2});
% termEC = MD_TERM(num_dims,{termEC_p,termEC_z});

terms = {termE1,termE2};
% terms = {termE1,termEA,termEB,termEC};


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


