function pde = advect_blob1(opts)
% 1D advect a blob across a domain. Use for feature tracking testing in the adaptiviy.  
%
% df/dt == -v * df/dx
%
% Run with
%
% asgard(@advect_blob1)
%
% Notes:
%  DLG - analytic solution won't work past about t = 0.002*500

%% Define the dimensions

dim_x = DIMENSION(-1,+1);
dimensions = {dim_x};
num_dims = numel(dimensions);

%% Define some parameters and add to pde object.

params.v = 2; % blob speed
params.sig = 0.1; % blob width

%% Define the analytic solution (optional).

a_x = @(x,p,t) exp(-(x-p.v*t).^2/p.sig.^2)+exp(-(x-p.v*t+2).^2/p.sig.^2);
a_t = @(t,p) t.*0+1;
sol1 = new_md_func(num_dims,{a_x,a_t});
solutions = {sol1};

%% Define the initial conditions

initial_conditions = {sol1};

%% Define the terms of the PDE
 
% -2 * df/dx

g1 = @(x,p,t,dat) x.*0-p.v;
pterm1 = GRAD(num_dims,g1,-1,'P','P');

term1_x = SD_TERM({pterm1});
term1   = MD_TERM(num_dims,{term1_x});

%%
% Add terms to the pde object

terms = {term1};

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


