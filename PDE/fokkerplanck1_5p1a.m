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
% explicit
% fk6d(fokkerplanck1_5p1a,5,3,0.01)
%
% implicit
% fk6d(fokkerplanck1_5p1a,5,4,3,[],[],1,'SG',[],1.5)

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
        a = 2;
        ret = 4.0/(sqrt(pi)*a^3) * exp(-x.^2/a^2);
    end

    function ret = x_psi(x) % manage the singularity at x=0      
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
dim_x.BCR = 'N'; % not set (equivalent to neumann)
dim_x.domainMin = 0;
dim_x.domainMax = +10;
dim_x.init_cond_fn = @(z,p) f0(z);

%%
% Add dimensions to the pde object
% Note that the order of the dimensions must be consistent with this across
% the remainder of this PDE.

pde.dimensions = {dim_x};

%% Setup the terms of the PDE

%% 
% d/dx*x^2*psi(x)/x*df/dx

term1_x.type = 'diff';
% Eq 1 : d/dx * x*psi(x) * q
term1_x.G1 = @(x,p,t,dat) x_psi(x);
term1_x.LF1 = -1; % Upwind
term1_x.BCL1 = 'D';
term1_x.BCR1 = 'D';
% Eq 2 : q = df/dx
term1_x.G2 = @(x,p,t,dat) x.*0+1;
term1_x.LF2 = +1; % Downwind
term1_x.BCL1 = 'N';
term1_x.BCR1 = 'N';

term1 = term_fill({term1_x});

%%
% d/dx*x^2*2*psi(x)*f

term2_x.type = 'grad';
term2_x.G1 = @(x,p,t,dat) x.^2*2.*psi(x);
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
    function dt=set_dt(pde)
        
        dims = pde.dimensions;
        xRange = dims{1}.domainMax-dims{1}.domainMin;
        lev = dims{1}.lev;
        CFL = pde.CFL;
        dx = xRange/2^lev;
        dt = CFL * dx;
        
    end

pde.set_dt = @set_dt;

end


