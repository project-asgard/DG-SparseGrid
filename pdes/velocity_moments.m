function pde_system = velocity_moments( opts )

%% Parameters:

n     = 1.0;
u     = 0.0;
theta = 0.4;

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
%% f = M(v):

%% Define the analytic solution (optional).

soln_x = @(x,p,t) 0*x+1;
soln_v = @(v,p,t) n/sqrt(2*pi*theta)*exp(-(v-u).^2./(2*theta));
soln_t = @(t,p)   0*t+1;
soln   = new_md_func(num_dims_xv,{soln_x,soln_v,soln_t});

%% Initial conditions

initial_conditions = {soln};

%% Analytic solution

analytic_solution = {soln};

%% Solution Vector

f = UNKNOWN( opts, dimensions_xv, analytic_solution, initial_conditions );

%
%% M = M[rho](v) (same as f):

%% Solution Vector

M = UNKNOWN( opts, dimensions_xv, analytic_solution, initial_conditions );

%
%% m_f_0 and m_M_0:

%% Define the analytic solution (optional).

soln_x = @(x,p,t) 0.*x+n;
soln_t = @(t,p)   0.*t+1;
soln   = new_md_func(num_dims_x,{soln_x,soln_t});

%% Initial conditions

initial_conditions = {soln};

%% Analytic solution

analytic_solution = {soln};

%% Solution Vector

m_f_0 = UNKNOWN( opts, dimensions_x, analytic_solution, initial_conditions );
m_M_0 = UNKNOWN( opts, dimensions_x, analytic_solution, initial_conditions );

%
%% m_f_1 and m_M_1:

%% Define the analytic solution (optional).

soln_x = @(x,p,t) 0.*x+n*u;
soln_t = @(t,p)   0.*t+1;
soln   = new_md_func(num_dims_x,{soln_x,soln_t});

%% Initial conditions

initial_conditions = {soln};

%% Analytic solution

analytic_solution = {soln};

%% Solution Vector

m_f_1 = UNKNOWN( opts, dimensions_x, analytic_solution, initial_conditions );
m_M_1 = UNKNOWN( opts, dimensions_x, analytic_solution, initial_conditions );

%
%% m_f_2 and m_M_2:

%% Define the analytic solution (optional).

soln_x = @(x,p,t) 0.*x+0.5*n*(u^2+theta);
soln_t = @(t,p)   0.*t+1;
soln   = new_md_func(num_dims_x,{soln_x,soln_t});

%% Initial conditions

initial_conditions = {soln};

%% Analytic solution

analytic_solution = {soln};

%% Solution Vector

m_f_2 = UNKNOWN( opts, dimensions_x, analytic_solution, initial_conditions );
m_M_2 = UNKNOWN( opts, dimensions_x, analytic_solution, initial_conditions );

%
%% f equation (dummy):

equation_f = EQUATION( f, {}, 'evolution', '' );

%
%% M equation:

MicroMacro = MICRO_MACRO_1X1V( opts, dimensions_xv );

% M(rho):

term_M = @( opts, Q, t ) MicroMacro.evaluate_rhs_Maxwellian( opts, Q, t );

term_M = TERM( M, {m_M_0,m_M_1,m_M_2}, {term_M}, false );

equation_M = EQUATION( M, {term_M}, 'closure', '' );
% equation_M = EQUATION( M, {}, 'evolution', '' );

%
%% m_f and m_M equations:

dV  = @(x,p,t) 0.*x+1;
One = @(x,p,t,dat) 0.*x+1;
e_0 = @(x,p,t) 0.*x+1;
e_1 = @(x,p,t) x;
e_2 = @(x,p,t) .5.*x.^2;

%% m_f_0 (zeroth moment of f):

term_x = SD_TERM({MASS(One,'','',dV)});
term_v = SD_TERM({VELOCITY_MOMENT_TERM(e_0,dim_v.moment_dV)});
term   = MD_TERM(num_dims_xv,{term_x,term_v});

term_m_f_0 = TERM( m_f_0, {f}, {term}, true );

equation_m_f_0 = EQUATION( m_f_0, {term_m_f_0}, 'closure', '' );

%% m_f_1 (first moment of f):

term_x = SD_TERM({MASS(One,'','',dV)});
term_v = SD_TERM({VELOCITY_MOMENT_TERM(e_1,dim_v.moment_dV)});
term   = MD_TERM(num_dims_xv,{term_x,term_v});

term_m_f_1 = TERM( m_f_1, {f}, {term}, true );

equation_m_f_1 = EQUATION( m_f_1, {term_m_f_1}, 'closure', '' );

%% m_f_2 (second moment of f)

term_x = SD_TERM({MASS(One,'','',dV)});
term_v = SD_TERM({VELOCITY_MOMENT_TERM(e_2,dim_v.moment_dV)});
term   = MD_TERM(num_dims_xv,{term_x,term_v});

term_m_f_2 = TERM( m_f_2, {f}, {term}, true );

equation_m_f_2 = EQUATION( m_f_2, {term_m_f_2}, 'closure', '' );

%% m_M_0 (zeroth moment of M):

equation_m_M_0 = EQUATION( m_M_0, {}, 'evolution', '' );

%% m_M_1 (first moment of M):

equation_m_M_1 = EQUATION( m_M_1, {}, 'evolution', '' );

%% m_M_2 (second moment of M)

equation_m_M_2 = EQUATION( m_M_2, {}, 'evolution', '' );

    function [ dt ] = set_dt( pde_system, CFL )
        dt = 1.0;
    end

%% Create the PDE System:

pde_system = PDE_SYSTEM( opts, {equation_f,equation_m_f_0,equation_m_f_1,equation_m_f_2,equation_M,equation_m_M_0,equation_m_M_1,equation_m_M_2}, @set_dt );

end