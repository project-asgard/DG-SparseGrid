function pde = fokkerplanck1_pitch_C(opts)
% 1D pitch angle collisional term 
% df/dt == d/dz ( (1-z^2) df/dz ) 
% 
% Here we use LDG for this second order system. We impose homogeneous
% Neumann BCs on f
%
% df/dt == d/dz( (1-z^2) df/dz ) becomes
%
% eq1 : df/dt ==  div( q(z) )   [div ,g(z)=1,BCL=D,BCR=D]
% eq2 :  q(z) == grad( f(z) )   [grad,g(z)=1,BCL=N,BCR=N]
%
% Here div(F\hat(z)) = d/dz(sqrt(1-z^2)F) and 
%            grad(F) = sqrt(1-z^2)df/dz\hat(z)
%
% Run with
%
% explicit
% asgard(@fokkerplanck1_pitch_C);
%
% implicit use in FP paper
% asgard(@fokkerplanck1_pitch_C,'timestep_method','matrix_exponential','dt',1,'lev',3,'deg',3);

    function ret = soln(z,t)
        
        h = [3,0.5,1,0.7,3,0,3];
        
        ret = zeros(size(z));
        
        P_0 = lin_legendre2(z,numel(h),2); % get matlab normalized legendres
        
        for l=1:numel(h)
            
            L = l-1;
            P_m = legendre(L,z); % Make sure Lin's legendre gives the matlab result.
            P2 = P_m(1,:)';
            P = P_0(:,l);
            assert(norm(P-P2)<1e-8); 
            
            ret = ret + h(l) * P * exp(-L*(L+1)*t);
            
        end
        
    end

%% Define the dimensions

dim_z = DIMENSION(-1,+1);
dim_z.moment_dV = @(x,p,t) 0*x+1;
dimensions = {dim_z};
num_dims = numel(dimensions);

%% Define the analytic solution (optional)

soln_z = @(z,p,t) soln(z,t);
soln1 = new_md_func(num_dims,{soln_z});
solutions = {soln1};

%% Define initial conditions

ic_z = @(z,p,t) soln(z,0);
ic1 = new_md_func(num_dims,{ic_z});
initial_conditions = {ic1};

%% Define the terms of the PDE
%
% Here we have 1 term1, having only nDims=1 (x) operators.

%% 
% termC
%
% d/dz( (1-z^2) df/dz )

dV_z = @(z,p,t,dat) sqrt(1-z.^2);
g1 = @(z,p,t,dat) 0*z+1;

pterm1  =  DIV(num_dims,g1,'',-1,'D','D','','','',dV_z);
pterm2  = GRAD(num_dims,g1,'',+1,'N','N','','','',dV_z);
termC_z = SD_TERM({pterm1,pterm2});
termC   = MD_TERM(num_dims,{termC_z});

terms = {termC};

%% Define some parameters

params.parameter1 = 0;
params.parameter2 = 1;

%% Define sources

sources = {};

%% Define function to set time step

    function dt=set_dt(pde,CFL)
        dims = pde.dimensions;
        % for Diffusion equation: dt = C * dx^2
        lev = dims{1}.lev;
        xMax = dims{1}.max;
        xMin = dims{1}.min;
        xRange = xMax-xMin;
        dx = xRange/2^lev;
        dt = CFL*dx^2;
    end

%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,[],initial_conditions,solutions);

end


