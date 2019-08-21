function pde = fokkerplanck1_5p1a
% Test the momentum dynamics for the RE problem for E = R = 0
%
% x^2 * df/dt == d/dx * x^2 ( psi(x)/x * df/dx + 2*psi(x)*f )
%             == d/dx*x^2*psi(x)/x*df/dx  +  d/dx*x^2*2*psi(x)*f
%
% (note the x^2 moved to the left side)
%
%
% Run with
%
% asgard(fokkerplanck1_5p1a,'imp licit',true,'lev',3,'num_steps',50,'CFL',1.5)

pde.CFL = 0.01;

%% Setup the dimensions
% 
% Here we setup a 1D problem (x)

    function ret = psi(x,t)
        
        phi = erf(x);
        dphi_dx = 2./sqrt(pi) * exp(-x.^2);

        ret = 1./(2*x.^2) .* (phi - x.*dphi_dx);
    end

    function ret = f0(x)
        
        test = 2;
        
        ret = zeros(size(x));
        switch test
            case 1
                a = 2;
                ret = 4.0/(sqrt(pi)*a^3) * exp(-x.^2/a^2);
            case 2
                for i=1:numel(x)
                    if x(i) <= 5
                        ret(i) = 3/5^3;
                    else
                        ret(i) = 0;
                    end
                end
                
        end
    end

    function ret = x_psi(x) % manage the singularity at x=0
        ret = zeros(size(x));
        for i=1:numel(x)
            if abs(x(i))<1e-5
                ret(i) = x(i).*0;
            else
                ret(i) = x(i).*psi(x(i));
            end
        end
    end

    function ret = soln(x,t)
        ret = 4/sqrt(pi) * exp(-x.^2);
    end

dim_x.domainMin = 0;
dim_x.domainMax = +10;
dim_x.init_cond_fn = @(z,p,t) f0(z);

%%
% Add dimensions to the pde object
% Note that the order of the dimensions must be consistent with this across
% the remainder of this PDE.

pde.dimensions = {dim_x};
num_dims = numel(pde.dimensions);

%% Setup the terms of the PDE

%%
% x^2 * df/dt (LHS non-identity coeff requires special treatment)

g1 = @(x,p,t,dat) x.^2;
pterm1 = MASS(g1);

termLHS_x = TERM_1D({pterm1});
termLHS   = TERM_ND(num_dims,{termLHS_x});

pde.termsLHS = {termLHS};

%% 
% d/dx*x^2*psi(x)/x*df/dx, split as two first order as 
%
% df/dt == d/dx g1(x) q(x)   [grad,g1(x)=x^2*psi(x)/x, BCL=D, BCR=N]
%  q(x) == d/dx g2(x) f(x)   [grad,g2(x)=1, BCL=N, BCR=D]

g1      = @(x,p,t,dat) x_psi(x);
g2      = @(x,p,t,dat) x.*0+1;
pterm1  = GRAD(num_dims,g1,-1,'D','N');
pterm2  = GRAD(num_dims,g2,+1,'N','D');
term1_x = TERM_1D({pterm1,pterm2});
term1   = TERM_ND(num_dims,{term1_x});

%%
% d/dx*x^2*2*psi(x)*f

g1      = @(x,p,t,dat) x*2.*x_psi(x);
pterm1  = GRAD(num_dims,g1,-1,'N','D');
term2_x = TERM_1D({pterm1});
term2   = TERM_ND(num_dims,{term2_x});

%%
% Add terms to the pde object

pde.terms = {term1,term2};

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
        xRange = dims{1}.domainMax-dims{1}.domainMin;
        lev = dims{1}.lev;
        dx = xRange/2^lev;
        dt = CFL * dx;
        
    end

pde.set_dt = @set_dt;

end


