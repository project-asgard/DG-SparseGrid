function pde = fokkerplanck2_complete(opts)
% Combining momentum and pitch angle dynamics
%
% Full PDE from the 2D runaway electron paper 
%
% d/dt f(p,z) == -div(flux_C + flux_E + flux_R)
%
% where 
%
% flux_C is flux due to electron-ion collisions
% flux_E is the flux due to E accleration
% flux_R is the flux due to radiation damping
%
% -div(flux_C) == termC1 + termC2 + termC3
% 
% termC1 == 1/p^2*d/dp*p^2*Ca*df/dp
% termC2 == 1/p^2*d/dp*p^2*Cf*f
% termC3 == Cb(p)/p^4 * d/dz( (1-z^2) * df/dz )
%
% -div(flux_E) == termE1 + termE2
%
% termE1 == -E*z*f(z) * 1/p^2 (d/dp p^2 f(p))
% termE2 == -E*p*f(p) * d/dz (1-z^2) f(z)
%
% -div(flux_R) == termR1 + termR2
%
% termR1 == 1/p^2 d/dp p^2 gamma(p) p / tau f(p) * (1-z^2) * f(z)
% termR2 == -1/(tau*gam(p)) f(p) * d/dz z(1-z^2) f(z)
%
% Run with
%
% explicit
% asgard(@fokkerplanck2_complete,'CFL',0.01,'case',1)
%
% implicit
% asgard(@fokkerplanck2_complete,'timestep_method','CN','num_steps',20,'CFL',1.0,'deg',3,'lev',4,'case',1)
%
% with adaptivity
% asgard(@fokkerplanck2_complete,'timestep_method','CN','num_steps',20,'CFL',1.0,'deg',3,'lev',4, 'adapt', true,'case',1)
%
% NOTES
%
% For Lin's case the command that should reproduce her results
% asgard(@fokkerplanck2_complete,'timestep_method','BE','lev',6,'deg',5,'dt',0.01,'num_steps',3000,'time_independent_A',true,'case',1)
%
% David's adjustment can be run as follows ...
% asgard(@fokkerplanck2_complete,'timestep_method','BE','lev',3,'deg',6,'dt',1,'num_steps',50,'grid_type','FG','time_independent_A',true,'case',1)
% or, run with the high order time integrator ...
% asgard(f@okkerplanck2_complete,'timestep_method','ode15s','lev',3,'deg',5,'dt',50,'num_steps',1,'grid_type','SG','case',1)

%%
% Define a few relevant functions

nuEE = 1;
vT = 1;
test = opts.case_;
switch test
    case 1
        delta = 0.042;
        Z = 1;
        E = 0.0025;
        tau = 10^5;
    case 2
        delta = 0.042;
        Z = 1;
        E = 0.25;
        tau = 10^5;
    case 3 
        delta = 0.042;
        Z = 1;
        E = 0.0025;
        tau = 10^5;
    case 4
        delta = 0.3;
        Z = 5;
        E = 0.4;
        tau = 10^5;
end
gamma = @(p)sqrt(1+(delta*p).^2);
vx = @(p)1/vT*(p./gamma(p));

Ca = @(p)nuEE*vT^2*(psi(vx(p))./vx(p));

Cb = @(p)1/2*nuEE*vT^2*1./vx(p).*(Z+phi(vx(p))-psi(vx(p))+delta^4*vx(p).^2/2);

Cf = @(p)2*nuEE*vT*psi(vx(p));

    function ret = phi(x)
        ret = erf(x);
    end

    function ret = psi(x,t)     
        dphi_dx = 2./sqrt(pi) * exp(-x.^2);
        ret = 1./(2*x.^2) .* (phi(x) - x.*dphi_dx);   
        ix = find(abs(x)<1e-5); % catch singularity at boundary
        ret(ix) = 0;
    end


    function ret = f0_z(x)
        test = opts.case_;
        ret = zeros(size(x));
        switch test
            case 1
                ret = x.*0+1;
            case 2
                ret = x.*0+1;
            case 3
                h = [3,0.5,1,0.7,3,0,3];
                
                for l=1:numel(h)
                    
                    L = l-1;
                    P_m = legendre(L,x); % Use matlab rather than Lin's legendre.
                    P = P_m(1,:)';
                    
                    ret = ret + h(l) * P;
                    
                end
            case 4
                ret = x.*0 + 1;
        end
    end

    function ret = f0_p(x)
        test = opts.case_;       
        ret = zeros(size(x));
        switch test
            
            case 1
                for i=1:numel(x)
                    if x(i) <= 5
                        ret(i) = 3/(2*5^3);
                    else
                        ret(i) = 0;
                    end
                end
                
            case 2 
                a = 2;
                ret = 2/(sqrt(pi)*a^3) * exp(-x.^2/a^2);
            case 3
                ret = 2/(3*sqrt(pi)) * exp(-x.^2);
            case 4
                N = 1000;
                h = 20/N;
                Q = 0;
                Fun = @(p)exp(-2/delta^2*sqrt(1+delta^2*p.^2));
                for i = 1:N
                    x0 = (i-1)*h;
                    x1 = i*h;
                    [xi,w] = lgwt(20,x0,x1);
                    Q = Q+sum(w.*Fun(xi).*xi.^2);
                 end
                ret = exp(-2/delta^2*sqrt(1+delta^2*x.^2))/(2*Q);
        end
    end


    function ret = soln_z(x,~)
        ret = x.*0+1;
    end
    function ret = soln_p(x,t)
        ret = 2/sqrt(pi) * exp(-x.^2);
    end

%% Setup the dimensions 

dim_p = DIMENSION(0.1,+20);
dim_p.init_cond_fn = @(x,p,t) f0_p(x);

dim_z = DIMENSION(-1,+1);
dim_z.init_cond_fn = @(x,p,t) f0_z(x);

dimensions = {dim_p,dim_z};
num_dims = numel(dimensions);

%% Define the terms of the PDE

%% -div(flux_C) == termC1 + termC2 + termC3

%% 
% termC1 == 1/p^2*d/dp*p^2*Ca*df/dp
%
% becomes 
%
% termC1 == g1(p) q(p)        [mass, g1(p) = 1/p^2,  BC N/A]
%   q(p) == d/dp g2(p) r(p)   [grad, g2(p) = p^2*Ca, BCL=D,BCR=N]
%   r(p) == d/dp g3(p) f(p)   [grad, g3(p) = 1,      BCL=N,BCR=D]

g1 = @(x,p,t,dat) 1./x.^2;
g2 = @(x,p,t,dat) x.^2.*Ca(x);
g3 = @(x,p,t,dat) x.*0+1; 

pterm1  = MASS(g1);
pterm2  = GRAD(num_dims,g2,+1,'D','D');
pterm3  = GRAD(num_dims,g3,-1,'N','N');

term1_p = TERM_1D({pterm1,pterm2,pterm3});
termC1  = TERM_ND(num_dims,{term1_p,[]});

%%
% termC2 == 1/p^2*d/dp*p^2*Cf*f
%
% becomes
%
% termC2 == g1(p) q(p)       [mass, g1(p)=1/p^2,  BC N/A]
%   q(p) == d/dp g2(p) f(p)  [grad, g2(p)=p^2*Cf, BCL=N,BCR=D]

g1 = @(x,p,t,dat) 1./x.^2;
g2 = @(x,p,t,dat) x.^2.*Cf(x);

pterm1  = MASS(g1);
pterm2  = GRAD(num_dims,g2,-1,'N','N');

term2_p = TERM_1D({pterm1,pterm2});
termC2   = TERM_ND(num_dims,{term2_p,[]});

%%
% termC3 == Cb(p)/p^4 * d/dz( (1-z^2) * df/dz )
%
% becomes
%
% termC3 == q(p) r(z)
%   q(p) == g1(p)            [mass, g1(p) = Cb(p)/p^4, BC N/A]
%   r(z) == d/dz g2(z) s(z)  [grad, g2(z) = 1-z^2,     BCL=D,BCR=D]
%   s(z) == d/dz g3(z) f(z)  [grad, g3(z) = 1,         BCL=N,BCR=N]

g1 = @(x,p,t,dat) Cb(x)./x.^4;
pterm1  = MASS(g1);
term3_p = TERM_1D({pterm1});

g2 = @(x,p,t,dat) (1-x.^2);
g3 = @(x,p,t,dat) x.*0+1;
pterm1  = GRAD(num_dims,g2,+1,'D','D');
pterm2  = GRAD(num_dims,g3,-1,'N','N');

term3_z = TERM_1D({pterm1,pterm2});

termC3 = TERM_ND(num_dims,{term3_p,term3_z});

%% -div(flux_E) == termE1 + termE2

% termE1 == -E*z*f(z) * 1/p^2 (d/dp p^2 f(p))
%        == r(z) * q(p)
%   r(z) == g1(z) f(z)       [mass, g1(z) = -E*z,  BC N/A]
%   q(p) == g2(p) u(p)       [mass, g2(p) = 1/p^2, BC N/A]
%   u(p) == d/dp g3(p) f(p)  [grad, g3(p) = p^2,   BCL=N,BCR=D]

% g1 = @(x,p,t,dat) -E.*x;
% g2 = @(x,p,t,dat) 1./x.^2;
% g3 = @(x,p,t,dat) x.^2;

g1 = @(x,p,t,dat) -x;
g2 = @(x,p,t,dat) 1./x.^2;
g3 = @(x,p,t,dat) E*x.^2;

pterm1   = MASS(g1);
termE1_z = TERM_1D({pterm1});

pterm1 = MASS(g2); 
pterm2 = GRAD(num_dims,g3,0,'N','N');% Lin's Setting

termE1_p = TERM_1D({pterm1,pterm2});

termE1 = TERM_ND(num_dims,{termE1_p,termE1_z});

% termE2 == -E*p*f(p) * d/dz (1-z^2) f(z)
%        == q(p) * r(z)
%   q(p) == g1(p) f(p)       [mass, g1(p) = -E*p,  BC N/A]
%   r(z) == d/dz g2(z) f(z)  [grad, g2(z) = 1-z^2, BCL=N,BCR=N]

% g1 = @(x,p,t,dat) -E.*x;
% g2 = @(x,p,t,dat) 1-x.^2;
g1 = @(x,p,t,dat) -1./x;
g2 = @(x,p,t,dat) E*(1-x.^2);

pterm1   = MASS(g1);
termE2_p = TERM_1D({pterm1});

pterm1   = GRAD(num_dims,g2,+1,'N','N');% Lin's Setting

termE2_z = TERM_1D({pterm1});

termE2 = TERM_ND(num_dims,{termE2_p,termE2_z});

%% -div(flux_R) == termR1 + termR2

% termR1 == 1/p^2 d/dp p^2 gamma(p) p / tau f(p) * (1-z^2) * f(z)
%        == q(p) * r(z)
%   q(p) == g1(p) u(p)       [mass, g1(p) = 1/p^2,                BC N/A]
%   u(p) == d/dp g2(p) f(p)  [grad, g2(p) = p^3 * gamma(p) / tau, BCL=N,BCR=D]
%   r(z) == g3(z) f(z)       [mass, g3(z) = 1-z^2,                BC N/A]

g1 = @(x,p,t,dat) 1./x.^2;
g2 = @(x,p,t,dat) x.^3 .* gamma(x) ./ tau;
g3 = @(x,p,t,dat) 1-x.^2;

pterm1   = MASS(g1);% This is not needed - by Lin
pterm2   = GRAD(num_dims,g2,1,'N','N');% Lin's Setting

termR1_p = TERM_1D({pterm1,pterm2});

pterm1   = MASS(g3);
termR1_z = TERM_1D({pterm1});

termR1   = TERM_ND(num_dims,{termR1_p,termR1_z});

% termR2 == -1/(tau*gam(p)) f(p) * d/dz z(1-z^2) f(z)
%        == q(p) * r(z)
%   q(p) == g1(p) f(p)       [mass, g1(p) = -1/(tau*gamma(p)),    BC N/A]
%   r(z) == d/dz g2(z) f(z)  [grad, g2(z) = z(1-z^2),             BCL=N,BCR=N]

g1 = @(x,p,t,dat) -1./(tau.*gamma(x));
g2 = @(x,p,t,dat) x.*(1-x.^2);

pterm1   = MASS(g1);
termR2_p = TERM_1D({pterm1});

pterm1   = GRAD(num_dims,g2,0,'N','N');% Lin's Setting

termR2_z = TERM_1D({pterm1});

termR2 = TERM_ND(num_dims,{termR2_p, termR2_z});

terms = {termC1, termC2, termC3, termE1, termE2, termR1, termR2};

%% Define some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;


%% Define sources

sources = {};

%% Define the analytic solution (optional).
% This requires nDims+time function handles.

analytic_solutions_1D = { ...
    @(x,p,t) soln_p(x,t), ...
    @(x,p,t) soln_z(x,t), ...
    @(t,p) 1
    };

%% Define function to set time step
    function dt=set_dt(pde,CFL)      
        dims = pde.dimensions;
        xRange = dims{1}.max-dims{1}.min;
        lev = dims{1}.lev;
        dx = xRange/2^lev;
        dt = CFL * dx;       
    end

%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,analytic_solutions_1D);

end


