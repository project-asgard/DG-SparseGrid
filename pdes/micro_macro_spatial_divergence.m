function pde_system = micro_macro_spatial_divergence( opts )

% d_t rho = 0             , rho = (rho_0,rho_1,rho_2)^T
% d_t u = - d_x F(rho)    , u = (u_0,u_1,u_2)^T         , u(t=0) = 0
% d_t f = - d_x vM[rho](v), f=f(x,v)                    , f(t=0) = 0
% m = < e f >             , e = (1,v,v^2/2)^T

%% Define the dimensionality:

dim_x = DIMENSION(0,1,opts.lev(1));
dim_x.moment_dV = @(x,p,t,dat) 0.*x+1;
dimensions_x    = {dim_x};
num_dims_x      = numel(dimensions_x);

dim_v = DIMENSION(-6,+6,opts.lev(2));
dim_v.moment_dV = @(v,p,t,dat) 0.*v+1;
dimensions_xv   = {dim_x,dim_v};
num_dims_xv     = numel(dimensions_xv);

%
%% Unknowns:

%
%% rho_0, rho_1, rho_2:
%% Specify n=n(x), u(x)=1, and theta(x)=1 so that rho_0=rho_1=rho_2

%% Define the analytic solution (optional).

soln_x = @(x,p,t) 1+.1*sin(2*pi*x);
soln_t = @(t,p)   0.*t+1;
soln   = new_md_func(num_dims_x,{soln_x,soln_t});

%% Initial conditions

initial_conditions = {soln};

%% Analytic solution

analytic_solution = {soln};

%% Solution Vector

rho_0 = UNKNOWN( opts, dimensions_x, analytic_solution, initial_conditions );
rho_1 = UNKNOWN( opts, dimensions_x, analytic_solution, initial_conditions );
rho_2 = UNKNOWN( opts, dimensions_x, analytic_solution, initial_conditions );

%
%% u_0, u_1, u_2:

%% Define the analytic solution (optional).

soln_x = @(x,p,t) 0.*x+0;
soln_t = @(t,p)   0.*t+1;
soln   = new_md_func(num_dims_x,{soln_x,soln_t});

%% Initial conditions

initial_conditions = {soln};

%% Analytic solution

analytic_solution = {};

u_0 = UNKNOWN( opts, dimensions_x, analytic_solution, initial_conditions );
u_1 = UNKNOWN( opts, dimensions_x, analytic_solution, initial_conditions );
u_2 = UNKNOWN( opts, dimensions_x, analytic_solution, initial_conditions );

%
%% f:

%% Define the analytic solution (optional).

soln_x = @(x,p,t) 0.*x+0;
soln_v = @(v,p,t) 0.*v+0;
soln_t = @(t,p)   0*t+1;
soln   = new_md_func(num_dims_xv,{soln_x,soln_v,soln_t});

%% Initial conditions

initial_conditions = {soln};

%% Analytic solution

analytic_solution = {};

%% Solution Vector

f = UNKNOWN( opts, dimensions_xv, analytic_solution, initial_conditions );

%
%% m_0, m_1, m_2:

%% Define the analytic solution (optional).

soln_x = @(x,p,t) 0.*x+0;
soln_t = @(t,p)   0.*t+0;
soln   = new_md_func(num_dims_x,{soln_x,soln_t});

%% Initial conditions

initial_conditions = {soln};

%% Analytic solution

analytic_solution = {};

%% Solution Vector

m_0 = UNKNOWN( opts, dimensions_x, analytic_solution, initial_conditions );
m_1 = UNKNOWN( opts, dimensions_x, analytic_solution, initial_conditions );
m_2 = UNKNOWN( opts, dimensions_x, analytic_solution, initial_conditions );

%
%% rho_0, rho_1, rho_2 equations (dummy):

equation_rho_0 = EQUATION( rho_0, {}, 'evolution', '' );
equation_rho_1 = EQUATION( rho_1, {}, 'evolution', '' );
equation_rho_2 = EQUATION( rho_2, {}, 'evolution', '' );

%
%% u_0, u_1, u_2 equations:

Euler = EULER_1D( opts, dimensions_x, 'PERIODIC', 'KiU', 1, 0.0 );

% - F_0(rho)_x:

F_0 = @( opts, Q, t ) Euler.evaluate_rhs_1( opts, Q, t );

term_u_0 = TERM( u_0, {rho_0,rho_1,rho_2}, {F_0}, false );

equation_u_0 = EQUATION( u_0, {term_u_0}, 'evolution', '' );

% - F_1(rho)_x:

F_1 = @( opts, Q, t ) Euler.evaluate_rhs_2( opts, Q, t );

term_u_1 = TERM( u_1, {rho_0,rho_1,rho_2}, {F_1}, false );

equation_u_1 = EQUATION( u_1, {term_u_1}, 'evolution', '' );

% - F_2(rho)_x:

F_2 = @( opts, Q, t ) Euler.evaluate_rhs_3( opts, Q, t );

term_u_2 = TERM( u_2, {rho_0,rho_1,rho_2}, {F_2}, false );

equation_u_2 = EQUATION( u_2, {term_u_2}, 'evolution', '' );

%
%% f equation:

MicroMacro = MICRO_MACRO_1X1V( opts, dimensions_xv );

% -(vM(rho))_x:

term_f = @( opts, Q, t ) MicroMacro.evaluate_rhs_vDotGradMaxwellian( opts, Q, t );

term_f = TERM( f, {rho_0,rho_1,rho_2}, {term_f}, false );

equation_f = EQUATION( f, {term_f}, 'evolution', '' );

%
%% m_0, m_1, m_2 equations:

dV  = @(x,p,t) 0.*x+1;
One = @(x,p,t,dat) 0.*x+1;
e_0 = @(x,p,t) 0.*x+1;
e_1 = @(x,p,t) x;
e_2 = @(x,p,t) .5.*x.^2;

%% m_0 (zeroth moment of f):

term_x = SD_TERM({MASS(One,'','',dV)});
term_v = SD_TERM({VELOCITY_MOMENT_TERM(e_0,dim_v.moment_dV)});
term   = MD_TERM(num_dims_xv,{term_x,term_v});

term_m_0 = TERM( m_0, {f}, {term}, true );

equation_m_0 = EQUATION( m_0, {term_m_0}, 'closure', '' );

%% m_1 (first moment of f):

term_x = SD_TERM({MASS(One,'','',dV)});
term_v = SD_TERM({VELOCITY_MOMENT_TERM(e_1,dim_v.moment_dV)});
term   = MD_TERM(num_dims_xv,{term_x,term_v});

term_m_1 = TERM( m_1, {f}, {term}, true );

equation_m_1 = EQUATION( m_1, {term_m_1}, 'closure', '' );

%% m_2 (second moment of f)

term_x = SD_TERM({MASS(One,'','',dV)});
term_v = SD_TERM({VELOCITY_MOMENT_TERM(e_2,dim_v.moment_dV)});
term   = MD_TERM(num_dims_xv,{term_x,term_v});

term_m_2 = TERM( m_2, {f}, {term}, true );

equation_m_2 = EQUATION( m_2, {term_m_2}, 'closure', '' );

    function [ dt ] = set_dt( pde_system, CFL )
        dt = 1.0;
    end

%% Create the PDE System:

pde_system = PDE_SYSTEM( opts, {equation_rho_0,equation_rho_1,equation_rho_2,equation_u_0,equation_u_1,equation_u_2,equation_f,equation_m_0,equation_m_1,equation_m_2}, @set_dt );

end