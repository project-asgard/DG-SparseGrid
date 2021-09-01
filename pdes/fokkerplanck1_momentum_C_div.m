function pde = fokkerplanck1_momentum_C_div(opts)

% Test the momentum dynamics for the RE problem for E = R = 0
%
% df/dt == 1/p^2*d/dp*p^2 ( psi(p)/p*df/dp + 2*psi(p)*f )
%       == 1/p^2*d/dp*p^2 * (psi(p)/p*df/dp)  +  1/p^2*d/dp*p^2 * (2*psi(p)*f)
%          \                              /    \                         /
%           --------- term1 --------------      -------- term2 ----------
%
% p is a spherical coordinate (p,th,ph), so the div and grad look like 
%
% div[] = 1/p^2 * d/dp * p^2[], grad[] = d/dp[]
%
% and the volument_element dV = p^2
%
% 
% split into two div(flux) terms (term1 and term2)
%
% term1 is done using SLDG defining eta1(p)=sqrt(psi(p)/p)
%
% eq1 :  df/dt == div(eta(p) * q)        [pterm1: div (g(p)=eta(p),+1, BCL=?, BCR=?)]
% eq2 :      q == eta(p) * grad(f)       [pterm2: grad(g(p)=eta(p),-1, BCL=D, BCR=N)]
%
% coeff_mat = pterm1.mat * pterm2.mat
%
% term2 is just a div
%
% eq1 :  df/dt == div(2*psi(p) * f)       [pterm1: div(g(p)=2*psi(p),+1, BCL=?, BCR=?]
%
% coeff_mat = pterm1.mat1
%
% Run with
%
% asgard(@fokkerplanck1_momentum_C_div,'timestep_method','BE','lev',3,'deg',4,'num_steps',50,'dt',1.0,'case',1)
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

%% Define some parameters and add to pde object.

params.parameter1 = 0;

%% Define the dimensions

dV_p = @(x,p,t,d) x.^2;
dim_p = DIMENSION(0,+10);
dim_p.moment_dV = dV_p;
dimensions = {dim_p};
num_dims = numel(dimensions);

%% Define the analytic solution (optional).

soln_p = @(x,p,t) soln(x);
soln1 = new_md_func(num_dims,{soln_p});
solutions = {soln1};

%% Define initial conditions

ic_p = @(z,p,t) f0(z);
ic1 = new_md_func(num_dims,{ic_p});
initial_conditions = {soln1};

%% LHS terms (mass only)

LHS_terms = {};

%% RHS terms

% term1 is done using SLDG defining eta(p)=sqrt(psi(p)/p)
%
% eq1 :  df/dt == div(eta(p) * q)        [pterm1: div (g1(p)=eta(p),+1, BCL=D, BCR=N)]
% eq2 :      q == eta(p) * grad(f)       [pterm2: grad(g2(p)=eta(p),-1, BCL=N, BCR=D)]

eta = @(p) sqrt(psi(p)./p);
g1 = @(x,p,t,dat) eta(x);
g2 = @(x,p,t,dat) eta(x);

pterm1 = DIV (num_dims,g1,'',+1,'D','N','','','',dV_p);
pterm2 = GRAD(num_dims,g2,'',-1,'N','D','','','',dV_p);

term1_p = SD_TERM({pterm1,pterm2}); % order here is as in equation
term1   = MD_TERM(num_dims,{term1_p});

% term2 is just a div
%
% eq1 :  df/dt == div(2*psi(p) * f)       [pterm1: div(g3(p)=2*psi(p),+1, BCL=N, BCR=D]

g3 = @(x,p,t,dat) 2 .* psi(x);

pterm3 = DIV(num_dims,g3,'',-1,'N','D','','','',dV_p);

term2_p = SD_TERM({pterm3});
term2   = MD_TERM(num_dims,{term2_p});

terms = {term1,term2};

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

pde = PDE(opts,dimensions,terms,LHS_terms,sources,params,@set_dt,[],initial_conditions,solutions);

end


