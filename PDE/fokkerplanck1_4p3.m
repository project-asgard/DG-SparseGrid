function pde = fokkerplanck1_4p3
% Problem 4.3 from the RE paper - radiation damping term  
% df/dt == -d/dz ( z(1-z^2)f )
%
% Run with
%
% explicit
% asgard(fokkerplanck1_4p3)
%
% implicit
% asgard(fokkerplanck1_4p3,'implicit',true)

pde.CFL = 0.01;
sig = 0.1;

%% Setup the dimensions
% 
% Here we setup a 1D problem (x)

    function ret = phi(z,t)
        ret = z.*exp(-t) ./ sqrt(1+(exp(-2*t)-1).*(z.^2));
    end
    function ret = f0(z)
        
        caseNumber = 4;
        shift = 0.36;
        
        switch caseNumber
            case 1
                f = exp(-z.^2/sig^2);
            case 2
                f = exp(-(z-shift).^2/sig^2);
            case 3
                f = exp(-(z+shift).^2/sig^2);
            case 4
                f = exp(-(z-shift).^2/sig^2) + exp(-(z+shift).^2/sig^2);
        end    
        ret = f;
    end
    function ret = soln(z,t)
        p = phi(z,t);
        t1 = p.*(1-p.^2);
        t2 = z.*(1-z.^2);
        t3 = f0(p);
        ret = t1./t2.*t3;
    end

dim_z.name = 'z';
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
% -d/dz ( z(1-z^2)f )

g1 = @(z,p,t,dat) -z.*(1-z.^2);
pterm1  = GRAD(num_dims,g1,-1,'D','D');
term1_z = TERM_1D({pterm1});
term1   = TERM_ND(num_dims,{term1_z});

%%
% Add terms to the pde object

pde.terms = {term1};

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
        Lmax = pde.dimensions{1}.domainMax;
        Lmin = pde.dimensions{1}.domainMin;
        LevX = pde.dimensions{1}.lev;
        dt = (Lmax-Lmin)/2^LevX*CFL;
    end

pde.set_dt = @set_dt;

end

