function pde = fokkerplanck2_6p1_withLHS(opts)
% Combining momentum and pitch angle dynamics
%
% Problems 6.1, 6.2, and 6.3 from the RE paper.
%
% p^2 * d/dt f(p,z,t) == termC1 + termC2 + termC3
%
% termC1 == d/dp*p^2*Ca*df/dp
% termC2 == d/dp*p^2*Cf*f
% termC2 == Cb(p)/p^2 * d/dz( (1-z^2) * df/dz )
%
% Run with
%
% explicit
% asgard(@fokkerplanck2_6p1_withLHS,'CFL',0.01,'case',1)
%
% implicit
% asgard(@fokkerplanck2_6p1_withLHS,'timestep_method','CN','num_steps',20,'CFL',1.0,'deg',3,'lev',4,'case',1)
%
% NOTES
% DLG - case 3 doesn't seem to be working, although it does work for the
% nonLHS version of this PDE?

%%
% Define a few relevant functions

nuEE = 1;
vT = 1;
delta = 0.042;
Z = 1;
E = 0.0025;
tau = 10^5;
gamma = @(p)sqrt(1+(delta*p).^2);
vx = @(p)1/vT*(p./gamma(p));
p_min = 0.005;

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
                
        end
    end


    function ret = soln_z(x,t)
        ret = x.*0+1;
    end
    function ret = soln_p(x,t)
        % From mathematica ...
        % a = 2;
        % Integrate[ 4/(Sqrt[Pi] a^3) Exp[-p^2/a^2], {p, p_min, 10}]
        switch p_min
            case 0
                Q = 0.5;        % for p_min = 0
            case 0.001
                Q = 0.499718;   % for p_min = 0.001
            case 0.005
                Q = 0.49859;    % for p_min = 0.005
            case 0.01              
                Q = 0.497179;   % for p_min = 0.01
            case 0.1
                Q = 0.471814;   % for p_min = 0.1
            case 0.2
                Q = 0.443769;   % for p_pim = 0.2
            case 0.5
                Q = 0.361837;   % for p_min = 0.5
            case 1.0
                Q = 0.23975;     % for p_min = 1.0
        end
        ret = Q * 4/sqrt(pi) * exp(-x.^2);
    end

%% Define the dimensions 

dim_p = DIMENSION(p_min,+10);
dim_p.init_cond_fn = @(x,p,t) f0_p(x);

dim_z = DIMENSION(-1,+1);
dim_z.init_cond_fn = @(x,p,t) f0_z(x);

dimensions = {dim_p,dim_z};
num_dims = numel(dimensions);

%% Define the terms of the PDE

%%
% LHS_term == p^2 d/dt f(p,z,t)

g1 = @(x,p,t,dat) x.^2;
pterm1 = MASS(g1);

LHS_term_p = TERM_1D({pterm1});
LHS_term   = TERM_ND(num_dims,{LHS_term_p,[]});

LHS_terms = {LHS_term};

%% 
% termC1 == d/dp*p^2*Ca*df/dp
%
% becomes 
%
% termC1 == d/dp g2(p) r(p)   [grad, g2(p) = p^2*Ca, BCL=D,BCR=N]        
%   r(p) == d/dp g3(p) f(p)   [grad, g3(p) = 1,      BCL=N,BCR=D]

g2 = @(x,p,t,dat) x.^2.*Ca(x);
g3 = @(x,p,t,dat) x.*0+1; 

pterm2  = GRAD(num_dims,g2,+1,'D','N');
pterm3  = GRAD(num_dims,g3,-1,'N','D');
term1_p = TERM_1D({pterm2,pterm3});
termC1  = TERM_ND(num_dims,{term1_p,[]});

%%
% termC2 == d/dp*p^2*Cf*f
%
% becomes
%
% termC2 == d/dp g2(p) f(p)  [grad, g2(p)=p^2*Cf, BCL=N,BCR=D]

g2 = @(x,p,t,dat) x.^2.*Cf(x);

pterm2  = GRAD(num_dims,g2,+1,'N','D');
term2_p = TERM_1D({pterm2});
termC2   = TERM_ND(num_dims,{term2_p,[]});

%%
% termC3 == Cb(p)/p^2 * d/dz( (1-z^2) * df/dz )
%
% becomes
%
% termC3 == q(p) r(z)
%   q(p) == g1(p)            [mass, g1(p) = Cb(p)/p^2, BC N/A]
%   r(z) == d/dz g2(z) s(z)  [grad, g2(z) = 1-z^2,     BCL=D,BCR=D]
%   s(z) == d/dz g3(z) f(z)  [grad, g3(z) = 1,         BCL=N,BCR=N]

g1 = @(x,p,t,dat) Cb(x)./x.^2;
pterm1  = MASS(g1);
term3_p = TERM_1D({pterm1});

g2 = @(x,p,t,dat) (1-x.^2);
g3 = @(x,p,t,dat) x.*0+1;
pterm1  = GRAD(num_dims,g2,+1,'D','D');
pterm2  = GRAD(num_dims,g3,-1,'N','N');
term3_z = TERM_1D({pterm1,pterm2});

termC3 = TERM_ND(num_dims,{term3_p,term3_z});

terms = {termC1,termC2,termC3};

%% Defin some parameters and add to pde object.
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

pde = PDE(opts,dimensions,terms,LHS_terms,sources,params,@set_dt,analytic_solutions_1D);

end


