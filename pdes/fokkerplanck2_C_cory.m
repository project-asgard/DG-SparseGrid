function pde = fokkerplanck2_C_cory(opts)
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
% asgard(@fokkerplanck2_C,'CFL',0.01,'case',1)
%
% implicit
% asgard(@fokkerplanck2_C,'timestep_method','BE','num_steps',20,'CFL',1.0,'deg',3,'lev',4,'case',1)
%
% case = 1 % step function in momentum
% case = 2 %
% case = 3 % 

params = fokkerplanck_parameters(opts);

%% Setup the dimensions

dim_x = DIMENSION(-1,1);
dim_x.moment_dV = @(x,p,t) 0*x+1;
dim_p = DIMENSION(0,+10);
dim_p.moment_dV = @(p,parm,t) p.^2;
dim_z = DIMENSION(-1,+1);
dim_z.moment_dV = @(z,p,t) 0*z+1;
dimensions = {dim_x,dim_p,dim_z};
num_dims = numel(dimensions);

%% Define the analytic solution (optional)

%soln_p = @(x,p,t) p.soln_p(x,p,t);
%soln_z = @(x,p,t) p.soln_z(x,p,t);
%soln1 = new_md_func(num_dims,{soln_p,soln_z});
solutions = {};

%% Define the initial conditions

ic_x = @(x,p,t) 0*x+1; %% To be specified
ic_p = @(x,p,t) 0*x+1;
ic_z = @(x,p,t) 0*x+1;
ic1 = new_md_func(num_dims,{ic_x,ic_p,ic_z});
initial_conditions = {ic1};

%% Define the boundary conditions

%% Define the terms of the PDE


%% div(flux_E) == termE1 + termE2

%% 

% termE1+E2 == div_x( pzf )
%
% DIV in p and MASS in z
%
% Since z changes sign on [-1,1] we have to split the term into two terms
% and upwind based on the sign of z.

% Surface jacobian for p constant is p^2
dV_x = @(x,p,t,dat) 0*x+1;
dV_p = @(x,p,t,dat) x.^2;
dV_z = @(x,p,t,dat) 0*x+1;

% Term for z>0, then |z| = z and we use standard upwinding.  

g1 = @(x,p,t,dat) 0*x-1;
pterm1   =  DIV(num_dims,g1,'',-1,'D','N','','','',dV_x); %flux = g*{f} + C*|g|*0.5[f] where C=-1.
term_adv_x = SD_TERM({pterm1});

% For hyperbolic problem: D - inflow, N - outflow/free boundary

g1 = @(x,p,t,dat) x;
pterm1   = MASS(g1,'','',dV_p);
term_adv_p = SD_TERM({pterm1});

g1 = @(x,p,t,dat) x.*(x > 0);
pterm1   = MASS(g1,'','',dV_z);
term_adv_z = SD_TERM({pterm1});

term_adv_pos = MD_TERM(num_dims,{term_adv_x,term_adv_p,term_adv_z});

%(g(x,p,z)u(x,p,z),v_x(x,p,z)) 
%g(x,p,z) = g(x)g(p)g(z)

% Term for z<0.  
% If z < 0, then |z| = -z.  
% Assume f(p,z) = f_p(p)f_z(z) and v(p,z) = v_p(p)v_z(z).
% Let <.,.>_p be the edge integral in p and (.,.)_z be the volume integral
% in z.  
% Standard upwinding jump flux becomes
%  <|Ez|[f],[v]>_{z,p} =  <|E|[f_p],[v_p]>_p(|z|f_z,v_z)_z
%                      = -<|E|[f_p],[v_p]>_p( z f_z,v_z)_z
% So we downwind in p when z is negative.

g1 = @(x,p,t,dat) 0*x-1;
pterm1   =  DIV(num_dims,g1,'',+1,'N','D','','','',dV_x); 
term_adv_x = SD_TERM({pterm1});

g1 = @(x,p,t,dat) x;
pterm1   = MASS(g1,'','',dV_p);
term_adv_p = SD_TERM({pterm1});

g1 = @(x,p,t,dat) x.*(x < 0);
pterm1   = MASS(g1,'','',dV_z);
term_adv_z = SD_TERM({pterm1});

term_adv_neg = MD_TERM(num_dims,{term_adv_x,term_adv_p,term_adv_z});


%%
% Pitch angle scattering
%
% LDG in q_z variable gives MASS in x and p, LDG in z
%
% becomes
% MASS in x,p  (first eqn)             [mass, g1 = 1, BC=NA] 
%              (secnd eqn)             [mass, g1 = 1, BC=NA] 
%
% LDG in z
%       div( A*q_z*\hat{z} )           [div , g1(p) = 1, BCL=D,BRC=N]
%   q(p) == A*grad(f)*\hat{z}          [grad, g2(p) = 1, BCL=N,BCR=D]
% A = 1

% Surface jacobian for z constant is dx*p*dp*sqrt(1-z^2)*dz
dV_x = @(x,p,t,dat) 0*x+1;
dV_p = @(x,p,t,dat) x;
dV_z = @(x,p,t,dat) sqrt(1-x.^2);

% MASS in x
g1 = @(x,p,t,dat) 0*x+1;
pterm1  = MASS(g1,'','',dV_x);
term_pitch_angle_scat_x = SD_TERM({pterm1,pterm1});

% MASS in p
g1 = @(x,p,t,dat) 0*x+1;
pterm1  = MASS(g1,'','',dV_p);
term_pitch_angle_scat_p = SD_TERM({pterm1,pterm1});

% LDG in z %Neumann boundary conditions
g1 = @(x,p,t,dat) x.*0+1;
pterm1  =  DIV(num_dims,g1,'',+1,'D','D','','','',dV_z);
pterm2  = GRAD(num_dims,g1,'',-1,'N','N','','','',dV_z);
term_pitch_angle_scat_z = SD_TERM({pterm1,pterm2});

term_pitch_angle_scat = MD_TERM(num_dims,{term_pitch_angle_scat_x,...
                                          term_pitch_angle_scat_p,...
                                          term_pitch_angle_scat_z});

                       
%% Finally collect all terms 

terms = {term_adv_pos,term_adv_neg,term_pitch_angle_scat};

%% Define sources

sources = {};

%% Define function to set time step
    function dt=set_dt(pde,CFL)
        dims = pde.dimensions;
        xRange = dims{3}.max-dims{3}.min;
        lev = dims{3}.lev;
        dx = xRange/2^lev;
        dt = CFL * dx;
    end

%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,[],initial_conditions,solutions);

end


