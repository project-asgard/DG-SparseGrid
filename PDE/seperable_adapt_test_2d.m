function pde = seperable_adapt_test_2d
% 2D test case using continuity equation, i.e.,
%
% df/dt == -df/dx
%
% Run with
%
% asgard(seperable_adapt_test_2d,'timestep_method','CN','adapt',true,'many_solution_capable',true)

pde.CFL = 0.1;

w = 1.0;
C = 1.0;

    function ret = offset(t)
        ret = sin(t);
    end
    function ret = expf(x,t)
        ret = exp( -( x - offset(t) ).^2 ./ w );
    end


%% Define the analytic solution (optional).
% This requires nDims+time function handles.

solution1 = { ...
    @(x,p,t)  x .* expf(x,t), ...
    @(y,p,t)  y.*0 + 1, ...
    @(t)      1 ...
    };

solution2 = { ...
    @(x,p,t)  expf(x,t), ...
    @(y,p,t)  -C.*sin(2*pi*y), ...
    @(t)      1 ...
    };

pde.solutions = {solution1,solution2};

%%
% Add initial conditions

pde.initial_conditions = {solution1,solution2};

%% Setup the dimensions
%
% Here we setup a 2D problem (x,y)

dim_x.domainMin = -10;
dim_x.domainMax = +10;

dim_y.domainMin = -1;
dim_y.domainMax = +1;

%%
% Add dimensions to the pde object
% Note that the order of the dimensions must be consistent with this across
% the remainder of this PDE.

pde.dimensions = {dim_x,dim_y};
num_dims = numel(pde.dimensions);

%% Setup the terms of the PDE
%
% Here we have 2 terms, having only nDims=2 (x,y) operators.

%%
% -df/dx which is 
%
% d/dx g1(x) f(x,y)          [grad,g1(x)=-1, BCL=P, BCR=P]

g1 = @(x,p,t,dat) x*0-1;
pterm1  = GRAD(num_dims,g1,0,'D','D');
term1_x = TERM_1D({pterm1});
term1   = TERM_ND(num_dims,{term1_x,[]});

%%
% Add terms to the pde object

pde.terms = {term1};

%% Construct some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

pde.params = params;

%% Add an arbitrary number of sources to the RHS of the PDE
% Each source term must have nDims + 1 (which here is 2+1 / (x,v) + time) functions describing the
% variation of each source term with each dimension and time.
% Here we define 3 source terms.

source1_t1p1 = { ...
    @(x,p,t) 2./w.*expf(x,t).*x.^2*cos(t), ... 
    @(y,p,t) y.*0+1, ...
    @(t)     1 ...
    };

source1_t1p2 = { ...
    @(x,p,t)  -2/w.*expf(x,t).*x*cos(t)*sin(t), ...
    @(y,p,t)  y.*0+1, ...
    @(t)      1 ...
    };

source1_t2p1 = { ...
    @(x,p,t)  expf(x,t), ...
    @(y,p,t)  y.*0+1, ...
    @(t)      1 ...
    };

source1_t2p2 = { ...
    @(x,p,t)  -2/w.*x.^2.*expf(x,t), ...
    @(y,p,t)  y.*0+1, ...
    @(t)      1 ...
    };

source1_t2p3 = { ...
    @(x,p,t)  2/w.*x.*expf(x,t).*sin(t), ...
    @(y,p,t)  y.*0+1, ...
    @(t)      1 ...
    };

source2_t1p1 = { ...
    @(x,p,t)  -2*C/w.*expf(x,t).*x.*cos(t), ...
    @(y,p,t)  sin(2*pi*y), ...
    @(t)      1 ...
    };

source2_t1p2 = { ...
    @(x,p,t)  2*C/w.*expf(x,t).*cos(t).*sin(t), ...
    @(y,p,t)  sin(2*pi*y), ...
    @(t)      1 ...
    };

source2_t2p1 = { ...
    @(x,p,t)  2*C/w.*expf(x,t).*x, ...
    @(y,p,t)  sin(2*pi*y), ...
    @(t)      1 ...
    };

source2_t2p2 = { ...
    @(x,p,t)  -2*C/w.*expf(x,t).*sin(t), ...
    @(y,p,t)  sin(2*pi*y), ...
    @(t)      1 ...
    };

%%
% Add sources to the pde data structure
pde.sources = {source1_t1p1,source1_t1p2,...
    source1_t2p1,source1_t2p2,source1_t2p3...
    source2_t1p1,source2_t1p2,...
    source2_t2p1,source2_t2p2};

%%
% function to set dt

    function dt=set_dt(pde,CFL)
        
        Lmax = pde.dimensions{1}.domainMax;
        Lmin = pde.dimensions{1}.domainMin;
        LevX = pde.dimensions{1}.lev;
        dt = (Lmax-Lmin)/2^LevX*CFL;
    end

pde.set_dt = @set_dt;

end


