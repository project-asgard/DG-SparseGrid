function pde_system = micro_macro_equilibrium_advection( opts )

%
%% Solving the micro-macro decomposition of f_t + (vf)_x = 0:
%  rho_t + F(rho)_x + < e q > = 0
%  g_t + q = - ( M(rho)_t + (vM(rho))_x )
%  q = (vg)_x

%
%% Macro Unknowns:

%% Define the dimensionality:

dim_x = DIMENSION(0,1,opts.lev(1));
dim_x.moment_dV = @(x,p,t,dat) 0*x+1;
dimensions_Macro = {dim_x};
num_dims_Macro = numel(dimensions_Macro);

%
%% First unknown (rho_0):

%% Define the analytic solution (optional).

soln_x = @(x,p,t) 1.+.5*sin(2*pi*(x-t));
soln_t = @(t,p)   0.*t+1;
soln   = new_md_func(num_dims_Macro,{soln_x,soln_t});

%% Initial conditions

initial_conditions = {soln};

%% Analytic solution

analytic_solution = {soln};

%% Solution Vector

rho_0 = UNKNOWN( opts, dimensions_Macro, analytic_solution, initial_conditions );

%
%% Second unknown (rho_1):

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
%% Third unknown (rho_2):

%% Define the analytic solution (optional).

soln_x = @(x,p,t) 0.5*(1.+.5*sin(2*pi.*(x-t))).*(1.+.1);
soln_t = @(t,p)   0.*t+1.;
soln   = new_md_func(num_dims_Macro,{soln_x,soln_t});

%% Initial conditions

initial_conditions = {soln};

%% Analytic solution

analytic_solution = {soln};

%% Solution Vector

rho_2 = UNKNOWN( opts, dimensions_Macro, analytic_solution, initial_conditions );

%
%% Micro Unknowns (g and q)

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

%% Unknown (q):

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

q = UNKNOWN( opts, dimensions_Micro, analytic_solution, initial_conditions );

%
%% Define the terms of the Macro system:

Euler = EULER_1D( opts, dimensions_Macro, 'PERIODIC', 'LLF', 1, 0.0 );

% F_1(rho)_x:

F_1 = @( opts, Q, t ) Euler.evaluate_rhs_1( opts, Q, t );

term_rho_0 = TERM( rho_0, {rho_0,rho_1,rho_2}, {F_1}, false );

% < e_0 q >_v:

dV = @(x,p,t,dat) 0*x+1;
e0 = @(x,p,t)     0*x+1;

term_q = VELOCITY_MOMENT_MD_TERM( e0, dV );

term_rho_0_q = TERM( rho_0, {q}, {term_q}, true );

equation_rho_0 = EQUATION( rho_0, {term_rho_0,term_rho_0_q}, 'evolution', '' );

% F_2(rho)_x:

F_2 = @( opts, Q, t ) Euler.evaluate_rhs_2( opts, Q, t );

term_rho_1 = TERM( rho_1, {rho_0,rho_1,rho_2}, {F_2}, false );

% < e_1 q >_v:

e_1 = @(x,p,t) x;

term_q = VELOCITY_MOMENT_MD_TERM( e_1, dV );

term_rho_1_q = TERM( rho_1, {q}, {term_q}, true );

equation_rho_1 = EQUATION( rho_1, {term_rho_1,term_rho_1_q}, 'evolution', '' );

% F_3(rho)_x:

F_3 = @( opts, Q, t ) Euler.evaluate_rhs_3( opts, Q, t );

term_rho_2 = TERM( rho_2, {rho_0,rho_1,rho_2}, {F_3}, false );

% < e_2 q >_v:

e_2 = @(x,p,t) .5*x.^2;

term_q = VELOCITY_MOMENT_MD_TERM( e_2, dV );

term_rho_2_q = TERM( rho_2, {q}, {term_q}, true );

equation_rho_2 = EQUATION( rho_2, {term_rho_2,term_rho_2_q}, 'evolution', '' );

%
%% Define the terms of the Micro system:

%% Terms for g equation.

MicroMacro = MICRO_MACRO_1X1V( opts, dimensions_Micro );

dV = @(x,p,t,dat) 0*x+1;
m1 = @(x,p,t,dat) 0*x-1;
p1 = @(x,p,t,dat) 0*x+1;

% q:

term_x = SD_TERM({MASS(p1,'','',dV)});
term_v = SD_TERM({MASS(p1,'','',dV)});
term_MM_1 = MD_TERM(num_dims_Micro,{term_x,term_v});

term_g_1 = TERM( g, {q}, {term_MM_1}, true );

% M(rho):

term_MM_2 = @( opts, Q, t ) MicroMacro.evaluate_rhs_Maxwellian( opts, Q, t );

term_g_2 = TERM( g, {rho_0,rho_1,rho_2}, {term_MM_2}, false );

% (vM(rho))_x:

term_MM_3 = @( opts, Q, t ) MicroMacro.evaluate_rhs_vDotGradMaxwellian( opts, Q, t );

term_g_3 = TERM( g, {rho_0,rho_1,rho_2}, {term_MM_3}, false );

% equation_g = EQUATION( g, {term_g_1,term_g_2,term_g_3}, 'evolution', '' );
equation_g = EQUATION( g, {term_g_1}, 'evolution', '' );%Temporarily turn off Maxwellian terms

%% Terms for q equation.

xp = @(x,p,t,dat) x.*(x>0);
xm = @(x,p,t,dat) x.*(x<0);

% - (vg)_x (for v>0):

term_x = SD_TERM({DIV(num_dims_Micro,m1,'',-1,'P','P','','','',dV)});
term_v = SD_TERM({MASS(xp,'','',dV)});
md_term_p = MD_TERM(num_dims_Micro,{term_x,term_v});

term_q_1 = TERM( q, {g}, {md_term_p}, true );

% - (vg)_x (for v<0):

term_x = SD_TERM({DIV(num_dims_Micro,m1,'',+1,'P','P','','','',dV)});
term_v = SD_TERM({MASS(xm,'','',dV)});
md_term_m = MD_TERM(num_dims_Micro,{term_x,term_v});

term_q_2 = TERM( q, {g}, {md_term_m}, true );

equation_q = EQUATION( q, {term_q_1,term_q_2}, 'closure', '' );

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

pde_system = PDE_SYSTEM( opts, {equation_rho_0,equation_rho_1,equation_rho_2,equation_g,equation_q}, @set_dt );

end