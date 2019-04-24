function pde = fokkerplanck1_4p4
% Problem 4.4 from the RE paper - evolution of the pitch angle dependence
% of f in the presence of electric field acceleration and collisions 
% 
% df/dt == -E d/dz((1-z^2) f) + C d/dz((1-z^2) df/dz)
%
% Run with
%
% explicit 
% fk6d(fokkerplanck1_4p4,4,2,0.2,[],[],0,[])
%
% implicit 
% fk6d(fokkerplanck1_4p4,4,4,1.5,[],[],1,[],[],0.2)

pde.CFL = 0.01;
sig = 0.1;
E = 4.0;
C = 1.0;

%% Setup the dimensions
% 
% Here we setup a 1D problem (x)

    function ret = phi(z,t)
        ret = z.*exp(-t) ./ sqrt(1-(exp(-2*t)-1).*(z.^2));
    end
    function ret = f0(z)
        
        caseNumber = 1;
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
        A = E/C;
        ret = A / (2*sinh(A)) * exp(A*z);
    end

BCL_fList = { ...
    @(z,p,t) 0, ...
    @(t,p) 0
    };

BCR_fList = { ...
    @(z,p,t) 0, ...
    @(t,p) 0
    };

dim_z.name = 'z';
dim_z.BCL = 'D'; % dirichlet
dim_z.BCL_fList = BCL_fList;
dim_z.BCR = 'D';
dim_z.BCR_fList = BCR_fList;
dim_z.domainMin = -1;
dim_z.domainMax = +1;
dim_z.lev = 2;
dim_z.deg = 2;
dim_z.FMWT = []; % Gets filled in later
dim_z.init_cond_fn = @(z,p) f0(z);

dim_z = checkDimension(dim_z);

%%
% Add dimensions to the pde object
% Note that the order of the dimensions must be consistent with this across
% the remainder of this PDE.

pde.dimensions = {dim_z};

%% Setup the terms of the PDE
%
% Here we have 1 term1, having only nDims=1 (x) operators.

%% 
% -E d/dz((1-z^2) f)

termE_z.type = 'grad'; % grad (see coeff_matrix.m for available types)
termE_z.G = @(z,p,t,dat) -E.*(1-z.^2); % G function for use in coeff_matrix construction.
termE_z.TD = 0; % Time dependent term or not.
termE_z.dat = []; % These are to be filled within the workflow for now
termE_z.LF = -1; % Upwind 
termE_z.name = 'd_dz';

termE = term_fill({termE_z});

%%
% +C d/dz((1-z^2) d/dz))

termC_z.type = 'diff';
% eq1 : 1 * dq/dx
termC_z.G1 = @(z,p,t,dat) z.*0+C;
termC_z.LF1 = -1; % upwind left
termC_z.BCL1 = 'N';
termC_z.BCR1 = 'N';
% eq2 : (1-z^2) * df/dx 
termC_z.G2 = @(z,p,t,dat) (1-z.^2);
termC_z.LF2 = +1; % upwind right
termC_z.BCL2 = 'D';
termC_z.BCR2 = 'D';

termC = term_fill({termC_z});


%%
% Add terms to the pde object

pde.terms = {termE,termC};

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

%% Other workflow options that should perhpas not be in the PDE?

pde.set_dt = @set_dt; % Function which accepts the pde (after being updated with CMD args).
pde.solvePoisson = 0; % Controls the "workflow" ... something we still don't know how to do generally. 
pde.applySpecifiedE = 0; % Controls the "workflow" ... something we still don't know how to do generally. 
pde.implicit = 0; % Can likely be removed and be a runtime argument. 
pde.checkAnalytic = 1; % Will only work if an analytic solution is provided within the PDE.

end

%% Define the various input functions specified above. 

%%
% Function to set time step
function dt=set_dt(pde)

Lmax = pde.dimensions{1}.domainMax;
LevX = pde.dimensions{1}.lev;
CFL = pde.CFL;
dt = Lmax/2^LevX*CFL;
end
