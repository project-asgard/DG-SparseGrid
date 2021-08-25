function pde = fokkerplanck2_E(opts)
% Combining momentum and pitch angle dynamics for the E term
%
% d/dt f(p,z) == div( F_p(f)\hat{p} + F_z(f)\hat{z} )
%
% where 
%
% F_p(f) = -Ezf    and    F_z(f) = -Esqrt(1-z^2)f
%
% Here
%
% div( A_p\hat{p} + A_z\hat{z} ) = 1/p^2*d/dp( p^2*A_p ) + 1/p*d/dz( sqrt(1-z^2)A_z )
%
% Run with
%
% asgard(@fokkerplanck2_E,'timestep_method','BE','dt',0.01,'num_steps',10,'case',3,'lev',4,'deg',4)
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
dim_p.moment_dV = @(x,p,t) x.^2;
dim_z.moment_dV = @(x,p,t) x.*0+1;
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

%% Boundary Conditions

switch opts.case_
    case 1
        BCL = soln1;
        BCR = soln1;
    case 2
        BCL = '';
        BCR = '';
    case 3
        BCL = soln1;
        BCR = soln1;
end

%% Define the terms of the PDE

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
pterm1   =  DIV(num_dims,g1,'',-1,'D','N',BCR,'','',dV_p);
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
pterm1   =  DIV(num_dims,g1,'',+1,'N','D','',BCL,'',dV_p);
termE2_p = SD_TERM({pterm1});

g1 = @(x,p,t,dat) x.*(x < 0);
pterm1   = MASS(g1,'','',dV_z);
termE2_z = SD_TERM({pterm1});

termE2 = MD_TERM(num_dims,{termE2_p,termE2_z});


%%%%%%
g1 = @(x,p,t,dat) 0*x-p.E;
pterm1   =  DIV(num_dims,g1,'',-1,'N','N','','','',dV_p);
termE4_p = SD_TERM({pterm1});

g1 = @(x,p,t,dat) x;
pterm1   = MASS(g1,'','',dV_z);
termE4_z = SD_TERM({pterm1});
termE4 = MD_TERM(num_dims,{termE4_p,termE4_z});

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
pterm1   =   DIV(num_dims,g1,'',-1,'D','D','','','',dV_z);
termE3_z = SD_TERM({pterm1});
termE3 = MD_TERM(num_dims,{termE3_p,termE3_z});

%%

terms = {termE1,termE2,termE3};
%terms = {termE4,termE3};

%% Define sources

switch opts.case_
    case 1
        sources = {};
    case 2
        sources = {};
    case 3
        ep2t = @(x,t) exp(-x.^2-t);
        
        fac=1;
        
        s1_p = @(x,p,t) -ep2t(x,t)*fac;
        s1_z = @(z,p,t) 1+cos(pi*z);
        s1 = new_md_func(num_dims,{s1_p,s1_z});
        
        s21_p = @(x,p,t) 2*ep2t(x,t)./x.*p.E*fac;
        s21_z = @(z,p,t) z;
        s21 = new_md_func(num_dims,{s21_p,s21_z});
        
        s22_p = @(x,p,t) -2*ep2t(x,t) .* x .* p.E*fac;
        s22_z = @(z,p,t) z;
        s22 = new_md_func(num_dims,{s22_p,s22_z});
        
        s23_p = @(x,p,t) +2*ep2t(x,t) ./ x .* p.E*fac;
        s23_z = @(z,p,t) z .* cos(pi*z);
        s23 = new_md_func(num_dims,{s23_p,s23_z});
        
        s24_p = @(x,p,t) -2*ep2t(x,t) .* x .* p.E*fac;
        s24_z = @(z,p,t) z .* cos(pi*z);
        s24 = new_md_func(num_dims,{s24_p,s24_z});
                
        s31_p = @(x,p,t) -2*ep2t(x,t) ./ x*fac .* p.E;
        s31_z = @(z,p,t) z;
        s31 = new_md_func(num_dims,{s31_p,s31_z});
        
        s32_p = @(x,p,t) -2*ep2t(x,t) ./ x*fac .* p.E;
        s32_z = @(z,p,t) z .* cos(pi*z);
        s32 = new_md_func(num_dims,{s32_p,s32_z});
        
        s33_p = @(x,p,t) -ep2t(x,t) ./ x*fac .* p.E;
        s33_z = @(z,p,t) pi * sin(pi * z);
        s33 = new_md_func(num_dims,{s33_p,s33_z});
        
        s34_p = @(x,p,t) ep2t(x,t) ./ x*fac .* p.E;
        s34_z = @(z,p,t) pi * z.^2 .* sin(pi * z);
        s34 = new_md_func(num_dims,{s34_p,s34_z});
        
        sources = {s1,s21,s22,s23,s24,s31,s32,s33,s34};
end

%% Construct moments

mass_moment_func = new_md_func(num_dims,{@(x,p,t) 0*x+1,...
                                    @(x,p,t) 0*x+1,...
                                    @(p,t)   0*t+1});
mass_moment = MOMENT({mass_moment_func});
moments = {mass_moment};

%% Define function to set time step
    function dt=set_dt(pde,CFL)
        dims = pde.dimensions;
        xRange = dims{1}.max-dims{1}.min;
        lev = dims{1}.lev;
        dx = xRange/2^lev;
        dt = CFL * dx;
    end

%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,[],initial_conditions,solutions,moments);

end


