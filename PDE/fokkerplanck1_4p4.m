function pde = fokkerplanck1_4p4
% Problem 4.4 from the RE paper - evolution of the pitch angle dependence
% of f in the presence of electric field acceleration and collisions 
% 
% df/dt == -E d/dz((1-z^2) f) + C d/dz((1-z^2) df/dz)
%
% Run with
%
% explicit 
% asgard(fokkerplanck1_4p4)
%
% implicit 
% asgard(fokkerplanck1_4p4,'implicit',true)

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
                f = exp(-z.^2/sig^2)/0.177245;
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

% dim_z.BCL = 'D'; % dirichlet
% dim_z.BCR = 'D';
dim_z.domainMin = -1;
dim_z.domainMax = +1;
dim_z.init_cond_fn = @(z,p,t) f0(z);

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
% -E d/dz((1-z^2) f)
% 
% termE_z.type = 'grad'; % grad (see coeff_matrix.m for available types)
% termE_z.G = @(z,p,t,dat) -E.*(1-z.^2); % G function for use in coeff_matrix construction.
% termE_z.LF = -1; % Upwind 
% 
% termE = term_fill({termE_z});

g1 = @(z,p,t,dat) -E.*(1-z.^2);
pterm1  = GRAD(num_dims,g1,-1,'D','D');
termE_z = TERM_1D({pterm1});
termE   = TERM_ND(num_dims,{termE_z});

%% 
% +C * d/dz( (1-z^2) df/dz )

% termC_z.type = 'diff';
% % eq1 : 1 * d/dx (1-z^2) q
% termC_z.G1 = @(z,p,t,dat) 1-z.^2;
% termC_z.LF1 = -1; % upwind left
% termC_z.BCL1 = 'D';
% termC_z.BCR1 = 'D';
% % eq2 : q = df/dx 
% termC_z.G2 = @(z,p,t,dat) z*0+1;
% termC_z.LF2 = +1; % upwind right
% termC_z.BCL2 = 'N';
% termC_z.BCR2 = 'N';
% 
% termC = term_fill({termC_z});

g1 = @(z,p,t,dat) 1-z.^2;
g2 = @(z,p,t,dat) z.*0+1;
pterm1  = GRAD(num_dims,g1,-1,'D','D');
pterm2  = GRAD(num_dims,g2,+1,'N','N');
termC_z = TERM_1D({pterm1,pterm2});
termC   = TERM_ND(num_dims,{termC_z});

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

