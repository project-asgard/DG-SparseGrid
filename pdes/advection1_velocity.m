function pde = advection1_velocity(opts)
%This is a test advection equation to help with the mirror velocity
%equation, which has the general form
%
% df/dt = 1/v^2 d/dv (v^3 f)
%
% Run with
%
%  asgard(@advection1_velocity,'timestep_method','BE', 'dt', 1, 'num_steps', 50, 'lev', 3, 'deg', 4)

% pde.CFL = 0.01;

%target properties
nu_s = 5; %slow down frequency in s^-1
m_a = 1.6726*10^-27; %mass of target in kg

%background properties
m_b = 9.109*10^-31; %mass of background in kg

soln = @(v,p,t) (m_a*nu_s/(m_a + m_b))*(v+t);

%% Define dimensions

dim_v = DIMENSION(0.01,10^2);
dim_v.moment_dV = @(x,p,t) x.^2;
dim_v.init_cond_fn = @(v,p,t) soln(v,p,t);

BCFunc_R = @(v,p,t) soln(v,p,t);
BCFunc_L = @(v,p,t) m_a*nu_s/(m_a + m_b);

% Domain is (a,b)

% The function is defined for the plane
% x = a and x = b
BCL_fList = { ...
    @(v,p,t) BCFunc_L(v,p,t), ... 
    @(t,p) t.*0 + 1
    };

BCR_fList = { ...
    @(v,p,t) BCFunc_R(v,p,t), ... % 
    @(t,p) t.*0 + 1
    };

%%
% Add dimensions to the pde object
% Note that the order of the dimensions must be consistent with this across
% the remainder of this PDE.

dimensions = {dim_v};
num_dims = numel(dimensions);

%% Define initial conditions

ic1 = new_md_func(num_dims,{soln});

initial_conditions = {ic1};

%% Define the terms of the PDE

% term V1 == 1/v^2 d/dv(v^3  f))
% term V1 == g(v) q(v)      [mass, g(v) = 1/v^2,  BC N/A]
% q(v) == d/dv(g2(v)f(v))   [grad, g2(v) = v^3, BCL= N, BCR=N]

dV = @(x,p,t,dat) x.^2;

g1 = @(v,p,t,dat) (m_a*nu_s/(m_a + m_b))*v;

pterm1  = DIV(num_dims,g1,'',-1,'N','D', '', BCR_fList,'',dV);
termV_s = SD_TERM({pterm1});
termV1   = MD_TERM(num_dims,{termV_s});

terms = {termV1};

%% Define some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;


%% Define sources
% Each source term must have nDims + 1

source = { ...
    @(v,p,t) (m_a.*nu_s/(m_a+m_b))*(1-m_a*nu_s/(m_a + m_b)*(4*v+3*t)),   ...   % s1v
    @(t,p) t.*0 + 1 ...   % s1t
    };

sources = {source};

%% Define the analytic solution (optional).

soln1 = new_md_func(num_dims,{soln});

solutions = {soln1};

%% Define to set time step
    function dt=set_dt(pde,CFL)      
        Lmax = pde.dimensions{1}.max;
        LevX = pde.dimensions{1}.lev;
        dt = Lmax/2^LevX*CFL;
    end

%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,[],initial_conditions,solutions);

end

