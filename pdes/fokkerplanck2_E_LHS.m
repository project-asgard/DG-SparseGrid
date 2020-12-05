function pde = fokkerplanck2_E_LHS(opts)
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
% asgard(@fokkerplanck2_E_LHS,'timestep_method','BE','dt',0.01,'num_steps',10,'case',3,'lev',4,'deg',4)
%
% case = 1 % flat initial condition, solution stays the same
% case = 2 % made up initial conditions, no known solution
% case = 3 % manufactured solution with sources

params = fokkerplanck_parameters(opts);

switch opts.case_
    case 3
        params.E = 0.25;
end

%% Setup the dimensions

dim_p = DIMENSION(0,+10);
dim_z = DIMENSION(-1,+1);
dim_p.jacobian = @(x,p,t) x.^2;
dim_z.jacobian = @(x,p,t) x.*0+1;
dimensions = {dim_p,dim_z};
num_dims = numel(dimensions);

%% Define the analytic solution (optional)

switch opts.case_
    case 1
        soln_p = @(x,p,t) x.*0+1;
        soln_z = @(x,p,t) x.*0+1;
        soln1 = new_md_func(num_dims,{soln_p,soln_z});
        solutions = {soln1};
    case 2
        solutions = {};
    case 3
        soln_p = @(x,p,t) exp(-x.^2);
        soln_z = @(z,p,t) cos(z.*pi)+1;
        soln_t = @(t,p) exp(-t);
        soln1 = new_md_func(num_dims,{soln_p,soln_z,soln_t});
        solutions = {soln1};
end


%% Define the initial conditions

switch opts.case_
    case 1
        ic_p = @(x,p,t) x.*0+1;
        ic_z = @(z,p,t) z.*0+1;
        ic1 = new_md_func(num_dims,{ic_p,ic_z});
    case 2
        %         ic_p = @(x,p,t) exp(-x.^2);
        %         ic_z = @(z,p,t) cos(z*pi)+1;
        ic_p = @(x,p,t) p.f0_p(x);
        ic_z = @(x,p,t) p.f0_z(x);
        ic1 = new_md_func(num_dims,{ic_p,ic_z});
    case 3
        ic1 = soln1;
end
initial_conditions = {ic1};


%% Define the terms of the PDE

% p^2 * df/dt (LHS non-identity coeff requires special treatment)

g1 = @(x,p,t,dat) x.^2;
pterm1 = MASS(g1);
LHS_term_p = SD_TERM({pterm1});
LHS_term   = MD_TERM(num_dims,{LHS_term_p,[]});
LHS_terms = {LHS_term};

%% -div(flux_E) == termE1 + termE2

% termE1 == -E*z*f(z) * (d/dp p^2 f(p))
%        == r(z) * q(p)
%   r(z) == g1(z) f(z)       [mass, g1(z) = -E*z,  BC N/A]
%   q(p) == g2(p) u(p)       [mass, g2(p) = 1, BC N/A]
%   u(p) == d/dp g3(p) f(p)  [grad, g3(p) = p^2,   BCL=N,BCR=D]

g1 = @(z,p,t,dat) -p.E .* z;
g3 = @(x,p,t,dat) x.^2;
pterm1   = MASS(g1);
termE1_z = SD_TERM({pterm1});
pterm2   = GRAD(num_dims,g3,0,'N','N'); % central flux due to change of z sign within domain
termE1_p = SD_TERM({pterm2});
termE1 = MD_TERM(num_dims,{termE1_p,termE1_z});

% figure(50)
% p = 0:.001:10;
% plot(p,g1(p,params,[],[]));
% hold on
% plot(p,g2(p,params,[],[]));
% plot(p,g3(p,params,[],[]));
% hold off

% termE2 == -E*p*f(p) * d/dz (1-z^2) f(z)
%        == q(p) * r(z)
%   q(p) == g1(p) f(p)       [mass, g1(p) = -E*p,  BC N/A]
%   r(z) == d/dz g2(z) f(z)  [grad, g2(z) = 1-z^2, BCL=N,BCR=N]

g1 = @(x,p,t,dat) -p.E .* x;
g2 = @(z,p,t,dat) 1-z.^2;
pterm1   = MASS(g1);
termE2_p = SD_TERM({pterm1});
pterm1   = GRAD(num_dims,g2,0,'N','N'); % central flux due to change of z sign within domain
termE2_z = SD_TERM({pterm1});
termE2= MD_TERM(num_dims,{termE2_p,termE2_z});

figure(51)
p = 0:.001:10;
plot(p,g1(p,params,[],[]));
hold on
plot(p,g2(p,params,[],[]));
hold off

terms = {termE1,termE2};

%% Define sources

switch opts.case_
    case 1
        sources = {};
    case 2
        sources = {};
    case 3
        ep2t = @(x,t) exp(-x.^2-t);
        
        fac=1;
        
        s1_p = @(x,p,t) -ep2t(x,t)*fac .* x.^2;
        s1_z = @(z,p,t) 1+cos(pi*z);
        s1 = new_md_func(num_dims,{s1_p,s1_z});
        
        s21_p = @(x,p,t) 2*ep2t(x,t)./x.*p.E*fac .* x.^2;
        s21_z = @(z,p,t) z;
        s21 = new_md_func(num_dims,{s21_p,s21_z});
        
        s22_p = @(x,p,t) -2*ep2t(x,t) .* x .* p.E*fac .* x.^2;
        s22_z = @(z,p,t) z;
        s22 = new_md_func(num_dims,{s22_p,s22_z});
        
        s23_p = @(x,p,t) +2*ep2t(x,t) ./ x .* p.E*fac .* x.^2;
        s23_z = @(z,p,t) z .* cos(pi*z);
        s23 = new_md_func(num_dims,{s23_p,s23_z});
        
        s24_p = @(x,p,t) -2*ep2t(x,t) .* x .* p.E*fac .* x.^2;
        s24_z = @(z,p,t) z .* cos(pi*z);
        s24 = new_md_func(num_dims,{s24_p,s24_z});
                
        s31_p = @(x,p,t) -2*ep2t(x,t) ./ x*fac .* p.E .* x.^2;
        s31_z = @(z,p,t) z;
        s31 = new_md_func(num_dims,{s31_p,s31_z});
        
        s32_p = @(x,p,t) -2*ep2t(x,t) ./ x*fac .* p.E .* x.^2;
        s32_z = @(z,p,t) z .* cos(pi*z);
        s32 = new_md_func(num_dims,{s32_p,s32_z});
        
        s33_p = @(x,p,t) -ep2t(x,t) ./ x*fac .* p.E .* x.^2;
        s33_z = @(z,p,t) pi * sin(pi * z);
        s33 = new_md_func(num_dims,{s33_p,s33_z});
        
        s34_p = @(x,p,t) ep2t(x,t) ./ x*fac .* p.E .* x.^2;
        s34_z = @(z,p,t) pi * z.^2 .* sin(pi * z);
        s34 = new_md_func(num_dims,{s34_p,s34_z});
        
        sources = {s1,s21,s22,s23,s24,s31,s32,s33,s34};
end

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


