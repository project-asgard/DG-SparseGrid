function pde = fokkerplanck2_complete
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
% termC3 == termC3 == Cb(p)/p^4 * d/dz( (1-z^2) * df/dz )
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
% asgard(fokkerplanck2_complete)
%
% implicit
% asgard(fokkerplanck2_complete,'implicit',true,'num_steps',20,'CFL',1.0,'deg',3,'lev',4)
%
% with adaptivity
% asgard(fokkerplanck2_complete,'implicit',true,'num_steps',20,'CFL',1.0,'deg',3,'lev',4, 'adapt', true)

pde.CFL = 0.01;

%%
% Select 6.1, 6.2, 6.3, etc where it goes as 6.test
test = '6p1b'; 

%%
% Define a few relevant functions

nuEE = 1;
vT = 1;
switch test
    case '6p1a'
        delta = 0.042;
        Z = 1;
        E = 0.0025;
        tau = 10^5;
    case '6p1b'
        delta = 0.042;
        Z = 1;
        E = 0.25;
        tau = 10^5;
    case '6p1c' 
        delta = 0.042;
        Z = 1;
        E = 0.0025;
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
        
        ret = zeros(size(x));
        switch test
            case '6p1a'
                ret = x.*0+1;
            case '6p1b'
                ret = x.*0+1;
            case '6p1c'
                h = [3,0.5,1,0.7,3,0,3];
                
                for l=1:numel(h)
                    
                    L = l-1;
                    P_m = legendre(L,x); % Use matlab rather than Lin's legendre.
                    P = P_m(1,:)';
                    
                    ret = ret + h(l) * P;
                    
                end
        end
    end

    function ret = f0_p(x)
               
        ret = zeros(size(x));
        switch test
            
            case '6p1a'
                for i=1:numel(x)
                    if x(i) <= 5
                        ret(i) = 3/(2*5^3);
                    else
                        ret(i) = 0;
                    end
                end
                
            case '6p1b' 
                a = 2;
                ret = 2/(sqrt(pi)*a^3) * exp(-x.^2/a^2);
            case '6p1c'
                ret = 2/(3*sqrt(pi)) * exp(-x.^2);
                
        end
    end


    function ret = soln_z(x,t)
        ret = x.*0+1;
    end
    function ret = soln_p(x,t)
        ret = 2/sqrt(pi) * exp(-x.^2);
    end

%% Setup the dimensions 

dim_p.domainMin = 0.1;
dim_p.domainMax = +7;
dim_p.init_cond_fn = @(x,p,t) f0_p(x);

dim_z.domainMin = -1;
dim_z.domainMax = +1;
dim_z.init_cond_fn = @(x,p,t) f0_z(x);

%%
% Add dimensions to the pde object
% Note that the order of the dimensions must be consistent with this across
% the remainder of this PDE.

pde.dimensions = {dim_p,dim_z};
num_dims = numel(pde.dimensions);

%% Setup the terms of the PDE

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
pterm2  = GRAD(num_dims,g2,+1,'N','N');
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

g1 = @(x,p,t,dat) -E.*x;
g2 = @(x,p,t,dat) 1./x.^2;
g3 = @(x,p,t,dat) x.^2;

pterm1   = MASS(g1);
termE1_z = TERM_1D({pterm1});

pterm1 = MASS(g2); 
pterm2 = GRAD(num_dims,g3,1,'N','N');% Lin's Setting
termE1_p = TERM_1D({pterm1,pterm2});

termE1 = TERM_ND(num_dims,{termE1_p,termE1_z});

% termE2 == -E*p*f(p) * d/dz (1-z^2) f(z)
%        == q(p) * r(z)
%   q(p) == g1(p) f(p)       [mass, g1(p) = -E*p,  BC N/A]
%   r(z) == d/dz g2(z) f(z)  [grad, g2(z) = 1-z^2, BCL=N,BCR=N]

g1 = @(x,p,t,dat) -E.*x;
g2 = @(x,p,t,dat) 1-x.^2;

pterm1   = MASS(g1);
termE2_p = TERM_1D({pterm1});

pterm1   = GRAD(num_dims,g2,0,'N','N');% Lin's Setting
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


%%
% Add terms to the pde object

pde.terms = {termC1, termC2, termC3, termE1, termE2, termR1, termR2};

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
    @(x,p,t) soln_p(x,t), ...
    @(x,p,t) soln_z(x,t), ...
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


