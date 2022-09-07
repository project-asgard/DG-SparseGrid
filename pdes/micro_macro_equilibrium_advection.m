function pde_system = micro_macro_equilibrium_advection( opts )

n     = @(x) 1.+.5*sin(2*pi*x);
u     = 1.0;
theta = 1.0;
Vmax  = 12.0;

%
%% Solving the micro-macro decomposition of f_t + (vf)_x = 0:
%  rho_t + F(rho)_x - < e q > = 0
%  g_t - q = - ( M(rho)_t + (vM(rho))_x )
%  q = - (vg)_x
%  rho_g = < e g >

%% Define the dimensionality:

dim_x = DIMENSION(0,1,opts.lev(1));
dim_x.moment_dV = @(x,p,t,dat) 0*x+1;
dimensions_x = {dim_x};
num_dims_x = numel(dimensions_x);

dim_v = DIMENSION(-Vmax,+Vmax,opts.lev(2));
dim_v.moment_dV = @(v,p,t,dat) 0*v+1;
dimensions_xv = {dim_x,dim_v};
num_dims_xv = numel(dimensions_xv);

%
%% Macro Unknowns:

%
%% First unknown (rho_0):

%% Define the analytic solution (optional).

soln_x = @(x,p,t) n(x-t);
soln_t = @(t,p)   0.*t+1;
soln   = new_md_func(num_dims_x,{soln_x,soln_t});

%% Initial conditions

initial_conditions = {soln};

%% Analytic solution

analytic_solution = {soln};

%% Solution Vector

rho_0 = UNKNOWN( opts, dimensions_x, analytic_solution, initial_conditions );

%
%% Second unknown (rho_1):

%% Define the analytic solution (optional).

soln_x = @(x,p,t) n(x-t).*u;
soln_t = @(t,p)   0.*t+1;
soln   = new_md_func(num_dims_x,{soln_x,soln_t});

%% Initial conditions

initial_conditions = {soln};

%% Analytic solution

analytic_solution = {soln};

%% Solution Vector

rho_1 = UNKNOWN( opts, dimensions_x, analytic_solution, initial_conditions );

%
%% Third unknown (rho_2):

%% Define the analytic solution (optional).

soln_x = @(x,p,t) 0.5*n(x-t).*(u^2+theta);
soln_t = @(t,p)   0.*t+1.;
soln   = new_md_func(num_dims_x,{soln_x,soln_t});

%% Initial conditions

initial_conditions = {soln};

%% Analytic solution

analytic_solution = {soln};

%% Solution Vector

rho_2 = UNKNOWN( opts, dimensions_x, analytic_solution, initial_conditions );

%
%% rho_g_0, rho_g_1, rho_g_2:

%% Define the analytic solution (optional).

soln_x = @(x,p,t) 0.*x+0;
soln_t = @(t,p)   0.*t+0;
soln   = new_md_func(num_dims_x,{soln_x,soln_t});

%% Initial conditions

initial_conditions = {soln};

%% Analytic solution

analytic_solution = {soln};

%% Solution Vector

rho_g_0 = UNKNOWN( opts, dimensions_x, analytic_solution, initial_conditions );
rho_g_1 = UNKNOWN( opts, dimensions_x, analytic_solution, initial_conditions );
rho_g_2 = UNKNOWN( opts, dimensions_x, analytic_solution, initial_conditions );

%
%% Micro Unknowns (g and q)

%% Unknown (g):

%% Define the analytic solution (optional).

soln_x = @(x,p,t) 0*x;
soln_v = @(v,p,t) 0*v;
soln_t = @(t,p)   0*t;
soln   = new_md_func(num_dims_xv,{soln_x,soln_v,soln_t});

%% Initial conditions

initial_conditions = {soln};

%% Analytic solution

analytic_solution = {soln};

%% Solution Vector

g = UNKNOWN( opts, dimensions_xv, analytic_solution, initial_conditions );

%% Unknown (q):

%% Define the analytic solution (optional).

soln_x = @(x,p,t) 0*x;
soln_v = @(v,p,t) 0*v;
soln_t = @(t,p)   0*t;
soln   = new_md_func(num_dims_xv,{soln_x,soln_v,soln_t});

%% Initial conditions

initial_conditions = {soln};

%% Analytic solution

analytic_solution = {soln};

%% Solution Vector

q = UNKNOWN( opts, dimensions_xv, analytic_solution, initial_conditions );

%
%% Define the terms of the Macro system:

dV  = @(x,p,t,dat) 0.*x+1;
One = @(x,p,t,dat) 0.*x+1;
e_0 = @(x,p,t)     0.*x+1;
e_1 = @(x,p,t)     x;
e_2 = @(x,p,t)     .5*x.^2;

Euler = EULER_1D( opts, dimensions_x, 'PERIODIC', 'KiU', 1, 0.0 );

%
%% rho_0 equation:

% F_1(rho)_x:

F_1 = @( opts, Q, t ) Euler.evaluate_rhs_1( opts, Q, t );

term_rho_0 = TERM( rho_0, {rho_0,rho_1,rho_2}, {F_1}, false, true );

% < e_0 q >_v:

term_x = SD_TERM({MASS(One,'','',dV)});
term_v = SD_TERM({VELOCITY_MOMENT_TERM(e_0,dim_v.moment_dV)});
term_q = MD_TERM(num_dims_xv,{term_x,term_v});

term_rho_0_q = TERM( rho_0, {q}, {term_q}, true, false );

equation_rho_0 = EQUATION( rho_0, {term_rho_0,term_rho_0_q}, 'evolution', '' );

%
%% rho_1 equation:

% F_2(rho)_x:

F_2 = @( opts, Q, t ) Euler.evaluate_rhs_2( opts, Q, t );

term_rho_1 = TERM( rho_1, {rho_0,rho_1,rho_2}, {F_2}, false, true );

% < e_1 q >_v:

term_x = SD_TERM({MASS(One,'','',dV)});
term_v = SD_TERM({VELOCITY_MOMENT_TERM(e_1,dim_v.moment_dV)});
term_q = MD_TERM(num_dims_xv,{term_x,term_v});

term_rho_1_q = TERM( rho_1, {q}, {term_q}, true, false );

equation_rho_1 = EQUATION( rho_1, {term_rho_1,term_rho_1_q}, 'evolution', '' );

%
%% rho_2 equation:

% F_3(rho)_x:

F_3 = @( opts, Q, t ) Euler.evaluate_rhs_3( opts, Q, t );

term_rho_2 = TERM( rho_2, {rho_0,rho_1,rho_2}, {F_3}, false, true );

% < e_2 q >_v:

term_x = SD_TERM({MASS(One,'','',dV)});
term_v = SD_TERM({VELOCITY_MOMENT_TERM(e_2,dim_v.moment_dV)});
term_q = MD_TERM(num_dims_xv,{term_x,term_v});

term_rho_2_q = TERM( rho_2, {q}, {term_q}, true, false );

equation_rho_2 = EQUATION( rho_2, {term_rho_2,term_rho_2_q}, 'evolution', '' );

%
%% Define the terms of the Micro system:

%% Terms for g equation.

MicroMacro = MICRO_MACRO_1X1V( opts, dimensions_xv );

dV = @(x,p,t,dat) 0*x+1;
m1 = @(x,p,t,dat) 0*x-1;
p1 = @(x,p,t,dat) 0*x+1;

% q:

term_x = SD_TERM({MASS(p1,'','',dV)});
term_v = SD_TERM({MASS(p1,'','',dV)});
term   = MD_TERM(num_dims_xv,{term_x,term_v});

term_g_1 = TERM( g, {q}, {term}, true, false );

% (vM(rho))_x:

term = @( opts, Q, t ) MicroMacro.evaluate_rhs_vDotGradMaxwellian( opts, Q, t );

term_g_2 = TERM( g, {rho_0,rho_1,rho_2}, {term}, false, true );

equation_g = EQUATION( g, {term_g_1,term_g_2}, 'evolution', '' );

%% Terms for q equation.

xp = @(x,p,t,dat) x.*(x>0);
xm = @(x,p,t,dat) x.*(x<0);

% - (vg)_x (for v>0):

term_x = SD_TERM({DIV(num_dims_xv,m1,'',-1,'P','P','','','',dV)});
term_v = SD_TERM({MASS(xp,'','',dV)});
term_p = MD_TERM(num_dims_xv,{term_x,term_v});

term_q_1 = TERM( q, {g}, {term_p}, true, false );

% - (vg)_x (for v<0):

term_x = SD_TERM({DIV(num_dims_xv,m1,'',+1,'P','P','','','',dV)});
term_v = SD_TERM({MASS(xm,'','',dV)});
term_m = MD_TERM(num_dims_xv,{term_x,term_v});

term_q_2 = TERM( q, {g}, {term_m}, true, false );

equation_q = EQUATION( q, {term_q_1,term_q_2}, 'closure', '' );

%% Equation for zeroth moment of g:

term_x = SD_TERM({MASS(One,'','',dV)});
term_v = SD_TERM({VELOCITY_MOMENT_TERM(e_0,dim_v.moment_dV)});
term   = MD_TERM(num_dims_xv,{term_x,term_v});

term_rho_g_0 = TERM( rho_g_0, {g}, {term}, true, false );

equation_rho_g_0 = EQUATION( rho_g_0, {term_rho_g_0}, 'closure', '' );

%% Equation for first moment of g:

term_x = SD_TERM({MASS(One,'','',dV)});
term_v = SD_TERM({VELOCITY_MOMENT_TERM(e_1,dim_v.moment_dV)});
term   = MD_TERM(num_dims_xv,{term_x,term_v});

term_rho_g_1 = TERM( rho_g_1, {g}, {term}, true, false );

equation_rho_g_1 = EQUATION( rho_g_1, {term_rho_g_1}, 'closure', '' );

%% Equation for second moment of g:

term_x = SD_TERM({MASS(One,'','',dV)});
term_v = SD_TERM({VELOCITY_MOMENT_TERM(e_2,dim_v.moment_dV)});
term   = MD_TERM(num_dims_xv,{term_x,term_v});

term_rho_g_2 = TERM( rho_g_2, {g}, {term}, true, false );

equation_rho_g_2 = EQUATION( rho_g_2, {term_rho_g_2}, 'closure', '' );

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

pde_system = PDE_SYSTEM( opts, {equation_rho_0,equation_rho_1,equation_rho_2,equation_g,equation_q,equation_rho_g_0,equation_rho_g_1,equation_rho_g_2}, @set_dt );

end