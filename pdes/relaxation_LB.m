function pde_system = relaxation_LB( opts )

%% Solving the PDE system:
% d_t f   = nu d_v ( (v-u[rho]) f + sqrt(theta[rho]) q )
% q       = sqrt(theta[rho])d_v f
% d_t rho = 0


n     = 1.0;
u     = 1.0;
theta = 0.5;
Vmax  = 6.0;

alpha_diffusion_flux = 0.5;

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
%% Unknowns:

%
%% f:

%% Define the analytic solution (optional).

soln_x = @(x,p,t) 0*x+1;
soln_v = @(v,p,t) exp(-(v-u).^2./(2*theta)).*n./sqrt(2*pi*theta);
soln_t = @(t,p)   0*t+1;
soln   = new_md_func(num_dims_xv,{soln_x,soln_v,soln_t});

%% Initial conditions

initial_conditions = {soln};

%% Analytic solution

analytic_solution = {soln};

%% Solution Vector

f = UNKNOWN( opts, dimensions_xv, analytic_solution, initial_conditions );

%
%% q:

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
%% rho_f_0:

%% Define the analytic solution (optional).

soln_x = @(x,p,t) 0.*x+n;
soln_t = @(t,p)   0.*t+1;
soln   = new_md_func(num_dims_x,{soln_x,soln_t});

%% Initial conditions

initial_conditions = {soln};

%% Analytic solution

analytic_solution = {soln};

%% Solution Vector

rho_f_0 = UNKNOWN( opts, dimensions_x, analytic_solution, initial_conditions );

%
%% rho_f_1:

%% Define the analytic solution (optional).

soln_x = @(x,p,t) 0.*x+n*u;
soln_t = @(t,p)   0.*t+1;
soln   = new_md_func(num_dims_x,{soln_x,soln_t});

%% Initial conditions

initial_conditions = {soln};

%% Analytic solution

analytic_solution = {soln};

%% Solution Vector

rho_f_1 = UNKNOWN( opts, dimensions_x, analytic_solution, initial_conditions );

%
%% rho_f_2:

%% Define the analytic solution (optional).

soln_x = @(x,p,t) 0.*x+0.5*n*(u^2+theta);
soln_t = @(t,p)   0.*t+1;
soln   = new_md_func(num_dims_x,{soln_x,soln_t});

%% Initial conditions

initial_conditions = {soln};

%% Analytic solution

analytic_solution = {soln};

%% Solution Vector

rho_f_2 = UNKNOWN( opts, dimensions_x, analytic_solution, initial_conditions );

LB = LENARD_BERNSTEIN_1X1V( opts, dimensions_xv, 'PERIODIC', 'ZEROFLUX', 0.0, [], alpha_diffusion_flux );

%% Terms for f equation.

term = @( opts, Q, t ) LB.evaluate_rhs_collision_operator_LB( opts, Q, t );

term_f = TERM( f, {f,q,rho_f_0,rho_f_1,rho_f_2}, {term}, false, true );

equation_f = EQUATION( f, {term_f}, 'evolution', '' );

%% Terms for q equation.

term = @( opts, Q, t ) LB.evaluate_rhs_diffusion_LB( opts, Q, t );

term_q = TERM( q, {f,rho_f_0,rho_f_1,rho_f_2}, {term}, false, true );

equation_q = EQUATION( q, {term_q}, 'closure', '' );

%% Terms for rho equations (dummy).

equation_rho_f_0 = EQUATION( rho_f_0, {}, 'evolution', '' );
equation_rho_f_1 = EQUATION( rho_f_1, {}, 'evolution', '' );
equation_rho_f_2 = EQUATION( rho_f_2, {}, 'evolution', '' );

    function [ dt ] = set_dt( pde_system, CFL )
        dim_v = pde_system.unknowns{1}.dimensions{2};
        vMax  = dim_v.max;
        vMin  = dim_v.min;
        dv    = ( vMax - vMin ) / ( 2^dim_v.lev );
        dt    = CFL * dv^2;
    end

%% Create the PDE System:

pde_system = PDE_SYSTEM( opts, {equation_f,equation_q,equation_rho_f_0,equation_rho_f_1,equation_rho_f_2}, @set_dt );

end