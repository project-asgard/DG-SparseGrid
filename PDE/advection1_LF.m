function pde = advection1_LF(opts)
% 1D test case where the advection coefficient switches sign within the domain [-pi,+pi]. 
% (see mathematica worksheet advection1_LF.nb)
%
% df/dt == -x * df/dx
%
% f(t=0) = Exp[-(x - 1)^2/sig] + Exp[-(x + 1)^2/sig];
% sig = 0.1;
%
% Run with
%
% asgard(@advection1_LF,'dt',0.01,'num_steps',100,'timestep_method','matrix_exponential','lev',4,'deg',4)
%
% Notes
% DLG - this works fine with central(0) or left(-1) flux, but fails with
% right(+1). Why?


%% Define the dimensions

dim_x = DIMENSION(-pi,+pi);
dimensions = {dim_x};
num_dims = numel(dimensions);

%% Define the analytic solution (optional).
    function ans = solution(x,p,t) 
        ep = (-1+x.*exp(-t)).^2;
        em = (+1+x.*exp(-t)).^2;
        ans = exp( -t -10*em -10*ep ) .* ( exp(10*em) + exp(10*ep) );
    end

a_x = @solution;
a_t = @(t,p) 1;
sol1 = new_md_func(num_dims,{a_x,a_t});
solutions = {sol1};

%% Define the initial conditions

ic_x = @solution;
ic1 = new_md_func(num_dims,{ic_x});
initial_conditions = {ic1};

%% Define the terms of the PDE

%% 
% -df/dx

g1 = @(x,p,t,dat) -x;
pterm1 = GRAD(num_dims,g1,0,'N','N');

term1_x = SD_TERM({pterm1});
term1   = MD_TERM(num_dims,{term1_x});

%%
% Add terms to the pde object

terms = {term1};

%% Define some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

%% Define sources

sources = {};

%% Define function to set time step

    function dt=set_dt(pde,CFL)
        
        dim = pde.dimensions{1};
        lev = dim.lev;
        xMax = dim.max;
        xMin = dim.min;
        xRange = xMax-xMin;
        dx = xRange/(2^lev);
        dt = CFL*dx;
    end

%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,[],initial_conditions,solutions);

end


