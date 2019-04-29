function pde = fokkerplanck2_6
% Combining momentum and pitch angle dynamics
%
% Problems 6.1, 6.2, and 6.3 from the RE paper.
%
% Run with
%
% explicit
% fk6d(fokkerplanck2_6p1,4,2,1)
%
% implicit
% fk6d(fokkerplanck2_6p1,4,2,100,[],[],1,[],[],5.0)

pde.CFL = 0.01;

%%
% Select 6.1, 6.2, 6.3, etc where it goes as 6.test
test = 3; 

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
        ret = 2/sqrt(pi) * exp(-x.^2);
    end

%% Setup the dimensions 

dim_p.BCL = 'N';
dim_p.BCR = 'D';
dim_p.domainMin = 0.5;
dim_p.domainMax = +10;
dim_p.init_cond_fn = @(x,p) soln_p(x,0);

dim_z.BCL = 'N';
dim_z.BCR = 'N';
dim_z.domainMin = -1;
dim_z.domainMax = +1;
dim_z.init_cond_fn = @(x,p) soln_z(x,0);

%%
% Add dimensions to the pde object
% Note that the order of the dimensions must be consistent with this across
% the remainder of this PDE.

pde.dimensions = {dim_p,dim_z};

%% Setup the terms of the PDE

%%
% x^2 * df/dt (LHS non-identity coeff requires special treatment)

termLHS_p.type = 'mass';
termLHS_p.G = @(x,p,t,dat) x.^2;

termLHS = term_fill({termLHS_p,[]});

pde.termsLHS = {termLHS};

%% 
% d/dp*p^2*Ca*df/dp

term1_p.type = 'diff';
% Eq 1 : d/dp * p^2*Ca * q
term1_p.G1 = @(x,p,t,dat) x.^2.*Ca(x);
term1_p.LF1 = +1; % Upwind
term1_p.BCL1 = 'D';
term1_p.BCR1 = 'N';
% Eq 2 : q = df/dx
term1_p.G2 = @(x,p,t,dat) x.*0+1;
term1_p.LF2 = -1; % Downwind
term1_p.BCL2 = 'N';
term1_p.BCR2 = 'D';

term1 = term_fill({term1_p,[]});

%%
% d/dp*p^2*Cf*f

term2_p.type = 'grad';
term2_p.G = @(x,p,t,dat) x.^2.*Cf(x);
term2_p.LF = +1;

term2 = term_fill({term2_p,[]});

%%
% Cb(p)/p^2 * d/dz( (1-z^2) * df/dz )

term3_p.type = 'mass';
term3_p.G = @(x,p,t,dat) Cb(x)./x.^2;

term3_z.type = 'diff';
% Eq 1 : d/dz (1-z^2) * q
term3_z.G1 = @(x,p,t,dat) (1-x.^2);
term3_z.LF1 = +1; % Upwind
term3_z.BCL1 = 'D';
term3_z.BCR1 = 'D';
% Eq 2 : q = df/dx
term3_z.G2 = @(x,p,t,dat) x.*0+1;
term3_z.LF2 = -1; % Downwind
term3_z.BCL2 = 'N';
term3_z.BCR2 = 'N';

term3 = term_fill({term3_p,term3_z});

%%
% Add terms to the pde object

pde.terms = {term1,term2,term3};
% pde.terms = {term1,term3};

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


