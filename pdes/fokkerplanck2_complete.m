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
% flux_C == A^2*F_C1 + F_C2
% 
% where 
%
% A\hat{p} = sqrt(C_A(p))\hat{p} and A\hat{z} = 1/p*sqrt(C_B(p))\hat{z}
%
% and
% 
% \F_C1 = grad(f) and \F_C2 = C_F(p)f\hat{p}
%
%
% flux_E == F_Ep(f)\hat{p} + F_Ez(f)\hat{z}
%
% where 
%
% F_Ep(f) = -Ezf    and    F_Ez(f) = -Esqrt(1-z^2)f
%
%
% flux_R == F_Rp \hat{p} + F_Rz \hat{z}
%                            
% where
%
% F_Rp = gamma(p)p/tau*(1-z^2)f(z) 
% F_Rz = -1/(tau*gamma(p))z*sqrt(1-z^2) f(z)
%
% Here
%
% div( A_p\hat{p} + A_z\hat{z} ) = 1/p^2*d/dp( p^2*A_p ) + 1/p*d/dz( sqrt(1-z^2)f )
%
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
% asgard(@fokkerplanck2_complete_div,'timestep_method','BE','lev',6,'deg',5,'dt',0.01,'num_steps',3000,'time_independent_A',true,'case',1)
%
% David's adjustment can be run as follows ...
% asgard(@fokkerplanck2_complete_div,'timestep_method','BE','lev',3,'deg',6,'dt',1,'num_steps',50,'grid_type','FG','time_independent_A',true,'case',1)
% or, run with the high order time integrator ...
% asgard(f@okkerplanck2_complete_div,'timestep_method','ode15s','lev',3,'deg',5,'dt',50,'num_steps',1,'grid_type','SG','case',1)

params = fokkerplanck_parameters(opts);

%% Setup the dimensions 

p_max = 10;
if isfield(opts.cmd_args,'p_max')
    p_max = opts.cmd_args.p_max;
end

dim_p = DIMENSION(0,p_max);
dim_z = DIMENSION(-1,+1);
dim_p.moment_dV = @(x,p,t,dat) x.^2;
dim_z.moment_dV = @(x,p,t,dat) x.*0+1;
dimensions = {dim_p,dim_z};
num_dims = numel(dimensions);

%% Define the initial conditions
ic_p = @(x,p,t) p.f0_p(x);
ic_z = @(x,p,t) p.f0_z(x);
%ic_p = @(x,p,t) 0*x+1; ic_z = @(x,p,t) 0*x+1;
ic1 = new_md_func(num_dims,{ic_p,ic_z});
initial_conditions = {ic1};

%% Define the analytic solution (optional)

solutions = {ic1};

%% Define the terms of the PDE

%% -div(flux_C) == termC1 + termC2 + termC3

%% C terms

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

%% E terms

%% div(flux_E) == termE1 + termE2

%% 

% termE1+E2 == div( F_p(f)\hat{p} )
%
% DIV in p and MASS in z
%
% Since z changes sign on [-1,1] we have to split the term into two terms
% and upwind based on the sign of z.

% Surface jacobian for p constant is p^2
dV_p = @(x,p,t,dat) x.^2;
dV_z = @(x,p,t,dat) 0*x+1;

% Term for z>0, then |z| = z and we use standard upwinding.  

g1 = @(x,p,t,dat) 0*x-p.E;
pterm1   =  DIV(num_dims,g1,'',-1,'D','N','','','',dV_p);
termE1_p = SD_TERM({pterm1});

g1 = @(x,p,t,dat) x.*(x > 0);
pterm1   = MASS(g1,'','',dV_z);
termE1_z = SD_TERM({pterm1});

termE1 = MD_TERM(num_dims,{termE1_p,termE1_z});

% Term for z<0.  
% If z < 0, then |z| = -z.  
% Assume f(p,z) = f_p(p)f_z(z) and v(p,z) = v_p(p)v_z(z).
% Let <.,.>_p be the edge integral in p and (.,.)_z be the volume integral
% in z.  
% Standard upwinding jump flux becomes
%  <|Ez|[f],[v]>_{z,p} =  <|E|[f_p],[v_p]>_p(|z|f_z,v_z)_z
%                      = -<|E|[f_p],[v_p]>_p( z f_z,v_z)_z
% So we downwind in p when z is negative.

g1 = @(x,p,t,dat) 0*x-p.E;
pterm1   =  DIV(num_dims,g1,'',+1,'N','D','','','',dV_p);
termE2_p = SD_TERM({pterm1});

g1 = @(x,p,t,dat) x.*(x < 0);
pterm1   = MASS(g1,'','',dV_z);
termE2_z = SD_TERM({pterm1});

termE2 = MD_TERM(num_dims,{termE2_p,termE2_z});

%%

% termE3 == div( F_z(f)\hat{z} )
%
% MASS in p and DIV in z
%
%

% Surface jacobian for z constant is p*sqrt(1-z^2)
dV_p = @(x,p,t,dat) x;
dV_z = @(x,p,t,dat) sqrt(1-x.^2);

g1 = @(x,p,t,dat) 0*x+1;
pterm1   =  MASS(g1,'','',dV_p);
termE3_p = SD_TERM({pterm1});

g1 = @(x,p,t,dat) -p.E*sqrt(1-x.^2);
pterm1   =   DIV(num_dims,g1,'',-1,'N','N','','','',dV_z);
termE3_z = SD_TERM({pterm1});
termE3 = MD_TERM(num_dims,{termE3_p,termE3_z});

%% R terms

%% div(flux_R) == termR1 + termR2

% termR1 == div( F_Rp \hat{p} )
%
% DIV in p and MASS in z

% Surface jacobian for p constant is p^2
dV_p = @(x,p,t,dat) x.^2;
dV_z = @(x,p,t,dat) 0*x+1;

g1 = @(x,p,t,dat) x .* p.gamma(x) ./ p.tau;
pterm1   =  DIV(num_dims,g1,'',-1,'N','D','','','',dV_p);
termR1_p = SD_TERM({pterm1});

g1 = @(x,p,t,dat) 1-x.^2; %Always positive
pterm1   = MASS(g1,'','',dV_z);
termR1_z = SD_TERM({pterm1});

termR1   = MD_TERM(num_dims,{termR1_p,termR1_z});

% termR2 == div( F_Rz \hat{z} )
%       
% MASS in p and DIV in z

% Surface jacobian for z constant is p*sqrt(1-z^2)
dV_p = @(x,p,t,dat) x;
dV_z = @(x,p,t,dat) sqrt(1-x.^2);

g1 = @(x,p,t,dat) x./(p.tau.*p.gamma(x));
pterm1   = MASS(g1,'','',dV_p);
termR2_p = SD_TERM({pterm1});

g1 = @(x,p,t,dat) -x.*sqrt(1-x.^2);
pterm1   =  DIV(num_dims,g1,'',-1,'N','N','','','',dV_z);
termR2_z = SD_TERM({pterm1});

termR2 = MD_TERM(num_dims,{termR2_p, termR2_z});

%terms = {termC1, termC2, termC3, termE1, termE2, termE3, termR1, termR2};
%terms = {termC1, termC2, termC3, termR1, termR2};
terms = {termC1, termC2, termC3, termE1, termE2, termE3};
%terms = {termC1, termC2, termC3};


%% Define sources

sources = {};

    function res = my_alpha(x,p,t)
%         disp(num2str(p.alpha_z(x)));
        res = p.alpha_z(x);
    end

switch opts.case_
    case 5
        source1_p = @(x,p,t,d) p.f0_p(x);
        source1_z = @(x,p,t,d) my_alpha(x,p,t);
        source1 = new_md_func(num_dims,{source1_p,source1_z});
        sources = {source1};
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

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,[],initial_conditions,solutions);

end


