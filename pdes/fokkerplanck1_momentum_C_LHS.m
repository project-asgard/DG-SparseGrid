function pde = fokkerplanck1_momentum_C_LHS(opts)
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
% asgard(@fokkerplanck1_momentum_C_LHS,'timestep_method','BE','lev',4,'deg',4,'num_steps',50,'CFL',1.5,'case',1)

    function ret = psi(x,t)
        
        phi = erf(x);
        dphi_dx = 2./sqrt(pi) * exp(-x.^2);

        ret = 1./(2*x.^2) .* (phi - x.*dphi_dx);
    end

    function ret = f0(x)
        
        test = opts.case_;
        
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

    function ret = soln(x)
        ret = 4/sqrt(pi) * exp(-x.^2);
    end

%% Define the dimensions

dim_x = DIMENSION(0,+10);
dimensions = {dim_x};
num_dims = numel(dimensions);

%% Define the initial conditions

ic_x = @(z,p,t) f0(z);
ic1 = new_md_func(num_dims,{ic_x});
initial_conditions = {ic1};

%% Define the terms of the PDE

% x^2 * df/dt (LHS non-identity coeff requires special treatment)

g1 = @(x,p,t,dat) x.^2;
pterm1 = MASS(g1);

LHS_term_x = SD_TERM({pterm1});
LHS_term   = MD_TERM(num_dims,{LHS_term_x});

LHS_terms = {LHS_term};

%% 
% d/dx*x^2*psi(x)/x*df/dx, split as two first order as 
%
% df/dt == d/dx g1(x) q(x)   [grad,g1(x)=x^2*psi(x)/x, BCL=D, BCR=N]
%  q(x) == d/dx g2(x) f(x)   [grad,g2(x)=1, BCL=N, BCR=D]

g1      = @(x,p,t,dat) x_psi(x);
g2      = @(x,p,t,dat) x.*0+1;
pterm1  = GRAD(num_dims,g1,-1,'D','N');
pterm2  = GRAD(num_dims,g2,+1,'N','D');
term1_x = SD_TERM({pterm1,pterm2});
term1   = MD_TERM(num_dims,{term1_x});

%%
% d/dx*x^2*2*psi(x)*f

g1      = @(x,p,t,dat) x*2.*x_psi(x);
pterm1  = GRAD(num_dims,g1,-1,'N','D');
term2_x = SD_TERM({pterm1});
term2   = MD_TERM(num_dims,{term2_x});

terms = {term1,term2};

%% Define some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;


%% Define sources

sources = {};

%% Define the analytic solution (optional).

soln_x = @(x,p,t) soln(x);
soln1 = new_md_func(num_dims,{soln_x});
solutions = {soln1};

%% Define function to set time step
    function dt=set_dt(pde,CFL)       
        dims = pde.dimensions;
        xRange = dims{1}.max-dims{1}.min;
        lev = dims{1}.lev;
        dx = xRange/2^lev;
        dt = CFL * dx;     
    end

%% Construct PDE

pde = PDE(opts,dimensions,terms,LHS_terms,sources,params,@set_dt,[],initial_conditions,solutions);

end


