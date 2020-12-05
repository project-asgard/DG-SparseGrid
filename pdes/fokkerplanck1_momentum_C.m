function pde = fokkerplanck1_momentum_C(opts)

% Test the momentum dynamics for the RE problem for E = R = 0
%
% df/dt == 1/x^2 * d/dx * x^2 ( psi(x)/x * df/dx + 2*psi(x)*f )
%       == 1/x^2 * d/dx*x^2*psi(x)/x*df/dx  +  1/x^2 * d/dx*x^2*2*psi(x)*f
%          \                              /    \                         /
%           --------- term1 --------------      -------- term2 ----------
% split into 
%
% term1
%
% eq1 :  df/dt == g(x) q(x)        [mass,g(x)=1/x^2,        BC N/A for mass]
% eq2 :   q(x) == d/dx g(x) p(x)   [grad,g(x)=x^2*psi(x)/x, BCL=D, BCR=N]
% eq3 :   p(x) == d/dx f(x)        [grad,g(x)=1,            BCL=N, BCR=D]
%
% coeff_mat = mat1 * mat2 * mat3
%
% term2
%
% eq1 :  df/dt == g(x) q(x)        [mass,g(x)=1/x^2,        BC N/A for mass]
% eq2 :   q(x) == d/dx g(x) f(x)   [grad,g(x)=x^2*2*psi(x), BCL=N, BCR=D]
%
% coeff_mat = mat1 * mat2
%
% Run with
%
% asgard(@fokkerplanck1_momentum_C,'timestep_method','CN','lev',3,'deg',4,'num_steps',50,'CFL',1.5,'case',1)
%
% case = 1 % maxwellian initial condition
% case = 2 % step function initial condition

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

%% Define the analytic solution (optional).

soln_z = @(z,p,t) soln(z);
soln1 = new_md_func(num_dims,{soln_z});
solutions = {soln1};

%% Define initial conditions

ic_x = @(z,p,t) f0(z);
ic1 = new_md_func(num_dims,{ic_x});
initial_conditions = {ic1};

%% Setup the terms of the PDE

% term1
%
% 1/x^2 * d/dx*x^2*psi(x)/x*df/dx
%
% eq1 :  df/dt == g1(x) q(x)        [mass,g1(x)=1/x^2,        BC N/A for mass]
% eq2 :   q(x) == d/dx g2(x) p(x)   [grad,g2(x)=x^2*psi(x)/x, BCL=D, BCR=N]
% eq3 :   p(x) == d/dx g3(x) f(x)   [grad,g3(x)=1,            BCL=N, BCR=D]
%
% coeff_mat = mat1 * mat2 * mat3

g1 = @(x,p,t,dat) 1./(x.^2);
g2 = @(x,p,t,dat) x_psi(x);
g3 = @(x,p,t,dat) x.*0+1;

pterm1 = MASS(g1);
pterm2 = GRAD(num_dims,g2,-1,'D','N');
pterm3 = GRAD(num_dims,g3,+1,'N','D');

term1_x = SD_TERM({pterm1,pterm2,pterm3});
term1   = MD_TERM(num_dims,{term1_x});

% term2
%
% 1/x^2 * d/dx*x^2*2*psi(x)*f
%
% eq1 :  df/dt == g(x) q(x)        [mass,g(x)=1/x^2,        BC N/A for mass]
% eq2 :   q(x) == d/dx g(x) f(x)   [grad,g(x)=x^2*2*psi(x), BCL=N, BCR=D]
%
% coeff_mat = mat1 * mat2

g1 = @(x,p,t,dat) 1./(x.^2);
g2 = @(x,p,t,dat) 2.*x.*x_psi(x);

pterm1 = MASS(g1);
pterm2 = GRAD(num_dims,g2,-1,'N','D');

term2_x = SD_TERM({pterm1,pterm2});
term2   = MD_TERM(num_dims,{term2_x});

terms = {term1,term2};

%% Define some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

%% Define sources

sources = {};

%% Define function to set time step
    function dt=set_dt(pde,CFL)       
        dims = pde.dimensions;
        xRange = dims{1}.max-dims{1}.min;
        lev = dims{1}.lev;
        dx = xRange/2^lev;
        dt = CFL * dx;       
    end

%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,[],initial_conditions,solutions);

end


