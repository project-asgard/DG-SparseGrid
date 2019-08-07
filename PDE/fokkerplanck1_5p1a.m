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
% asgard(fokkerplanck1_5p1a,'implicit',true,'lev',3,'num_steps',50,'CFL',1.5)

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

dim_x.BCL = 'N'; % neumann
dim_x.BCR = 'D'; % not set (equivalent to neumann)
dim_x.domainMin = 0;
dim_x.domainMax = +10;
dim_x.init_cond_fn = @(z,p,t) f0(z);

%%
% Add dimensions to the pde object
% Note that the order of the dimensions must be consistent with this across
% the remainder of this PDE.

pde.dimensions = {dim_x};

%% Setup the terms of the PDE

%%
% x^2 * df/dt (LHS non-identity coeff requires special treatment)

termLHS_x.type = 'mass';
termLHS_x.G = @(x,p,t,dat) x.^2;

termLHS = term_fill({termLHS_x});

pde.termsLHS = {termLHS};


%% 
% d/dx*x^2*psi(x)/x*df/dx

term1_x.type = 'diff';
% Eq 1 : d/dx * x*psi(x) * q
term1_x.G1 = @(x,p,t,dat) x_psi(x);
term1_x.LF1 = -1; % Upwind
term1_x.BCL1 = 'D';
term1_x.BCR1 = 'N';
% Eq 2 : q = df/dx
term1_x.G2 = @(x,p,t,dat) x.*0+1;
term1_x.LF2 = +1; % Downwind
term1_x.BCL2 = 'N';
term1_x.BCR2 = 'D';

term1 = term_fill({term1_x});


%%
% d/dx*x^2*2*psi(x)*f

term2_x.type = 'grad';
term2_x.G = @(x,p,t,dat) x*2.*x_psi(x);
term2_x.LF = -1;

term2 = term_fill({term2_x});


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


