function pde = advection2_LF(opts)
% 1D test case where the advection coefficient switches sign within the domain [-pi,+pi]. 
% (see mathematica worksheet advection1_LF.nb)
%
% df/dt == -x * df/dx - y * df/dy
%
% and is simply the tensor product of the advection1_LF case because I
% couldn't figure out how to solve it for real in 2D
%
% Run with
%
% implicit
% asgard(@advection2_LF,'dt',1e-2,'num_steps',15,'timestep_method','BE','lev',4,'deg',4)
%
% explicit
% asgard(@advection2_LF,'dt',1e-2,'num_steps',15,'timestep_method','RK3','lev',4,'deg',4)
%
% Notes
% DLG - this works fine with central(0) or left(-1) flux, but fails with
% right(+1). Why?


%% Define the dimensions

dim_x = DIMENSION(-pi,+pi);
dim_y = DIMENSION(-pi,+pi);
dimensions = {dim_x,dim_y};
num_dims = numel(dimensions);

%% Define the analytic solution (optional).
    function ans = solution(x,p,t) 
        ep = (-1+x.*exp(-t)).^2;
        em = (+1+x.*exp(-t)).^2;
        ans = exp( -t -10*em -10*ep ) .* ( exp(10*em) + exp(10*ep) );
    end

a_x = @solution;
a_y = a_x;
a_t = @(t,p) 1;
sol1 = new_md_func(num_dims,{a_x,a_y,a_t});
solutions = {sol1};

%% Define the initial conditions

ic_x = @solution;
ic_y = ic_x;
ic1 = new_md_func(num_dims,{ic_x,ic_y});
initial_conditions = {ic1};

%% Define the terms of the PDE

%% 
% -x*df/dx

g1 = @(x,p,t,dat) -x;
pterm1 = GRAD(num_dims,g1,0,'N','N');

term1_x = SD_TERM({pterm1});
term1   = MD_TERM(num_dims,{term1_x,[]});

%% 
% -y*df/dy

g1 = @(y,p,t,dat) -y;
pterm1 = GRAD(num_dims,g1,0,'N','N');

term2_y = SD_TERM({pterm1});
term2   = MD_TERM(num_dims,{[],term2_y});

%%
% Add terms to the pde object

terms = {term1,term2};

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


