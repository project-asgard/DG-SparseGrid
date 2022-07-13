function pde_system = micro_macro_equilibrium_advection( opts )

%
%% Solving a reduced version of the micro-macro system:
%  rho_t + F(rho)_x = 0
%  g_t = - ( M(rho)_t + (vM(rho))_x )

%
%% Macro Unknowns:

%% Define the dimensionality:

dim_x = DIMENSION(0,1,opts.lev(1));
dim_x.moment_dV = @(x,p,t,dat) 0*x+1;
dimensions_Macro = {dim_x};
num_dims_Macro = numel(dimensions_Macro);

%
%% First unknown (rho_1):

%% Define the analytic solution (optional).

soln_x = @(x,p,t) 1.+.5*sin(2*pi*(x-t));
soln_t = @(t,p)   0.*t+1;
soln   = new_md_func(num_dims_Macro,{soln_x,soln_t});

%% Initial conditions

initial_conditions = {soln};

%% Analytic solution

analytic_solution = {soln};

%% Solution Vector

rho_1 = UNKNOWN( opts, dimensions_Macro, analytic_solution, initial_conditions );

%
%% Second unknown (rho_2):

%% Define the analytic solution (optional).

soln_x = @(x,p,t) 1.+.5*sin(2*pi*(x-t));
soln_t = @(t,p)   0.*t+1;
soln   = new_md_func(num_dims_Macro,{soln_x,soln_t});

%% Initial conditions

initial_conditions = {soln};

%% Analytic solution

analytic_solution = {soln};

%% Solution Vector

rho_2 = UNKNOWN( opts, dimensions_Macro, analytic_solution, initial_conditions );

%
%% Third unknown (rho_3):

%% Define the analytic solution (optional).

soln_x = @(x,p,t) 0.5*(1.+.5*sin(2*pi.*(x-t))).*(1.+.1);
soln_t = @(t,p)   0.*t+1.;
soln   = new_md_func(num_dims_Macro,{soln_x,soln_t});

%% Initial conditions

initial_conditions = {soln};

%% Analytic solution

analytic_solution = {soln};

%% Solution Vector

rho_3 = UNKNOWN( opts, dimensions_Macro, analytic_solution, initial_conditions );

%
%% Micro Unknown

%% Define the dimensionality:

dim_v = DIMENSION(-6,+6,opts.lev(2));
dim_v.moment_dV = @(v,p,t,dat) 0*v+1;
dimensions_Micro = {dim_x,dim_v};
num_dims_Micro = numel(dimensions_Micro);

%% Unknown (g):

%% Define the analytic solution (optional).

soln_x = @(x,p,t) 0*x;
soln_v = @(v,p,t) 0*v;
soln_t = @(t,p)   0*t;
soln   = new_md_func(num_dims_Micro,{soln_x,soln_v,soln_t});

%% Initial conditions

initial_conditions = {soln};

%% Analytic solution

analytic_solution = {soln};

%% Solution Vector

g = UNKNOWN( opts, dimensions_Micro, analytic_solution, initial_conditions );

%
%% Define the terms of the Macro system:

Euler = EULER_1D( opts, dimensions_Macro, 'PERIODIC', 'LLF', 1, 0.0 );

% F_1(rho)_x:

F_1 = @( opts, Q, t ) Euler.evaluate_rhs_1( opts, Q, t );

descriptor = {F_1};

term_rho_1 = TERM( rho_1, {rho_1,rho_2,rho_3}, descriptor, false );

equation_rho_1 = EQUATION( rho_1, {term_rho_1}, 'evolution', '' );

% F_2(rho)_x:

F_2 = @( opts, Q, t ) Euler.evaluate_rhs_2( opts, Q, t );

descriptor = {F_2};

term_rho_2 = TERM( rho_2, {rho_1,rho_2,rho_3}, descriptor, false );

equation_rho_2 = EQUATION( rho_2, {term_rho_2}, 'evolution', '' );

% F_3(rho)_x:

F_3 = @( opts, Q, t ) Euler.evaluate_rhs_3( opts, Q, t );

descriptor = {F_3};

term_rho_3 = TERM( rho_3, {rho_1,rho_2,rho_3}, descriptor, false );

equation_rho_3 = EQUATION( rho_3, {term_rho_3}, 'evolution', '' );

%
% Define the terms of the Macro system:

MicroMacro = MICRO_MACRO_1X1V( opts, dimensions_Micro );

% M(rho):

term_MM_1 = @( opts, Q, t ) MicroMacro.evaluate_rhs_Maxwellian( opts, Q, t );

term_g_1 = TERM( g, {rho_1,rho_2,rho_3}, {term_MM_1}, false );

% (vM(rho))_x:

term_MM_2 = @( opts, Q, t ) MicroMacro.evaluate_rhs_vDotGradMaxwellian( opts, Q, t );

term_g_2 = TERM( g, {rho_1,rho_2,rho_3}, {term_MM_2}, false );

equation_g = EQUATION( g, {term_g_1,term_g_2}, 'evolution', '' );

    function [ dt ] = set_dt( pde_system, CFL )
        deg   = pde_system.unknowns{4}.deg;
        dim_x = pde_system.unknowns{4}.dimensions{1};
        dim_v = pde_system.unknowns{4}.dimensions{2};
        xMax  = dim_x.max;
        xMin  = dim_x.min;
        vMax  = dim_v.max;
        vMin  = dim_v.min;
        dx    = ( xMax - xMin ) / ( 2^dim_x.lev );
        dt    = CFL * dx / max([abs(vMax),abs(vMin)]) / ( 2 * deg - 1 );
    end

%% Create the PDE System:

pde_system = PDE_SYSTEM( opts, {equation_rho_1,equation_rho_2,equation_rho_3,equation_g}, @set_dt );

end