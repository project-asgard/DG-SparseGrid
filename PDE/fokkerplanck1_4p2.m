function pde = fokkerplanck1_4p2
% 1D pitch angle collisional term 
% df/dt == d/dz ( (1-z^2) df/dz ) 
% 
% Here we use LDG for this second order system. We impose homogeneous
% Neumann BCs on f
%
% df/dt == d/dz( (1-z^2) df/dz ) becomes
%
% eq1 : df/dt == d/dz g(z) q(z)       [grad,g(z)=1-z^2,BCL=D,BCR=D]
% eq2 :  q(z) == d/dz g(z) f(z)       [grad,g(z)=1,BCL=N,BCR=N]
%
%
% Run with
%
% explicit
% asgard(fokkerplanck1_4p2);
%
% implicit
% asgard(fokkerplanck1_4p2,'implicit',true,'num_steps',20,'lev',3,'deg',3);

pde.CFL = 0.005;

%% Setup the dimensions
% 
% Here we setup a 1D problem (x)

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

dim_z.domainMin = -1;
dim_z.domainMax = +1;
dim_z.init_cond_fn = @(z,p,t) soln(z,0);

%%
% Add dimensions to the pde object
% Note that the order of the dimensions must be consistent with this across
% the remainder of this PDE.

pde.dimensions = {dim_z};
num_dims = numel(pde.dimensions);

%% Setup the terms of the PDE
%
% Here we have 1 term1, having only nDims=1 (x) operators.

%% 
% termC
%
% d/dz( (1-z^2) df/dz )

g1 = @(z,p,t,dat) 1-z.^2;
g2 = @(z,p,t,dat) z.*0+1;
pterm1  = GRAD(num_dims,g1,-1,'D','D');
pterm2  = GRAD(num_dims,g2,+1,'N','N');
termC_z = TERM_1D({pterm1,pterm2});
termC   = TERM_ND(num_dims,{termC_z});

%%
% Add terms to the pde object

pde.terms = {termC};

%% Construct some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

pde.params = params;

%% Add an arbitrary number of sources to the RHS of the PDE
% Each source term must have nDims + 1

%%
% Add sources to the pde data structure
pde.sources = {};

%% Define the analytic solution (optional).
% This requires nDims+time function handles.

pde.analytic_solutions_1D = { ...
    @(z,p,t) soln(z,t), ...
    @(t,p) 1 
    };

%%
% Function to set time step

    function dt=set_dt(pde,CFL)
        dims = pde.dimensions;
        % for Diffusion equation: dt = C * dx^2
        lev = dims{1}.lev;
        xMax = dims{1}.domainMax;
        xMin = dims{1}.domainMin;
        xRange = xMax-xMin;
        dx = xRange/2^lev;
        dt = CFL*dx^2;
    end

pde.set_dt = @set_dt;

end

