function pde = fokkerplanck1_4p4(opts)
% Problem 4.4 from the RE paper - evolution of the pitch angle dependence
% of f in the presence of electric field acceleration and collisions 
% 
% df/dt == -E d/dz((1-z^2) f) + C d/dz((1-z^2) df/dz)
%
% Run with
%
% explicit 
% asgard(@fokkerplanck1_4p4,'CFL',0.01,'case',3)
%
% implicit 
% asgard(@fokkerplanck1_4p4,'timestep_method','BE','num_steps',50,'dt',0.05,'lev',4,'deg',4,'case',1)

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
        
        caseNumber = opts.case_;
        shift = 0.36;
        scaling = 0.177245;
        
        switch caseNumber
            case 1
                f = exp(-z.^2/sig^2)/scaling;
            case 2
                f = exp(-(z-shift).^2/sig^2)/scaling;
            case 3
                f = exp(-(z+shift).^2/sig^2)/scaling;
            case 4
                f = (exp(-(z-shift).^2/sig^2) + exp(-(z+shift).^2/sig^2))/(2*scaling);
        end    
        ret = f;
    end
    function ret = soln(z,t)
        A = E/C;
        ret = A / (2*sinh(A)) * exp(A*z);
    end


dim_z = DIMENSION(-1,+1);
dim_z.init_cond_fn = @(z,p,t) f0(z);

dimensions = {dim_z};
num_dims = numel(dimensions);

%% Define the terms of the PDE
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

terms = {termE,termC};

%% Define some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;


%% Define sources

sources = {};

%% Define the analytic solution (optional).
% This requires nDims+time function handles.

analytic_solutions_1D = { ...
    @(z,p,t) soln(z,t), ...
    @(t,p) 1 
    };

%% Define function to set time step

    function dt=set_dt(pde,CFL)
        
        Lmax = pde.dimensions{1}.max;
        Lmin = pde.dimensions{1}.min;
        LevX = pde.dimensions{1}.lev;
        dt = (Lmax-Lmin)/2^LevX*CFL;
    end

%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,analytic_solutions_1D);

end

