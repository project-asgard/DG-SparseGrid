function pde_system = euler_system1( opts )

switch opts.case_
    case 1
        TEST = 'LinearAdvection';
    case 2
        TEST = 'Riemann';
    otherwise
        TEST = 'LinearAdvection';
end

assert(any(strcmp(TEST,{'LinearAdvection','Riemann'})))

%
%% Solving the Euler equations as a system:
%  u_t + f(u)_x = 0

%% Define the dimensionality:

dim_x = DIMENSION(0,1,opts.lev(1));
dim_x.moment_dV = @(x,p,t,dat) 0*x+1;
dimensions = {dim_x};
num_dims = numel(dimensions);

%% First unknown (density: u_1):

%% Define the analytic solution (optional).

switch TEST
    case 'LinearAdvection'
        soln_x = @(x,p,t) 1.+.5*sin(2*pi*(x-t));
    case 'Riemann'
        soln_x = @(x,p,t) 0.125 + (1.0-0.125).*(x <= 0.5);
end
soln_t = @(t,p)   0.*t+1;
soln   = new_md_func(num_dims,{soln_x,soln_t});

%% Initial conditions

initial_conditions = {soln};

%% Analytic solution

analytic_solution = {soln};

%% Solution Vector

u_1 = UNKNOWN( opts, dimensions, analytic_solution, initial_conditions );

%% Second unknown (momentum: u_2):

%% Define the analytic solution (optional).
switch TEST
    case 'LinearAdvection'
        soln_x = @(x,p,t) 1.+.5*sin(2*pi*(x-t));
    case 'Riemann'
        soln_x = @(x,p,t) 0.*x;
end
soln_t = @(t,p)   0.*t+1;
soln   = new_md_func(num_dims,{soln_x,soln_t});

%% Initial conditions

initial_conditions = {soln};

%% Analytic solution

analytic_solution = {soln};

%% Solution Vector

u_2 = UNKNOWN( opts, dimensions, analytic_solution, initial_conditions );

%% Third unknown (energy: u_3):

%% Define the analytic solution (optional).
switch TEST
    case 'LinearAdvection'
        soln_x = @(x,p,t) 0.5*(1.+.5*sin(2*pi.*(x-t))).*(1.+1.e-10);
    case 'Riemann'
        soln_x = @(x,p,t) 0.05 + (0.5-0.05).*(x <= 0.5);
end
soln_t = @(t,p)   0.*t+1.;
soln   = new_md_func(num_dims,{soln_x,soln_t});

%% Initial conditions

initial_conditions = {soln};

%% Analytic solution

analytic_solution = {soln};

%% Solution Vector

u_3 = UNKNOWN( opts, dimensions, analytic_solution, initial_conditions );

%% Define the terms of the PDE system

switch TEST
    case 'LinearAdvection'
        BCs   = 'PERIODIC';
        T_min = 0.0;
    case 'Riemann'
        BCs   = 'HOMOGENEOUS';
        T_min = 1.0e-8;
end

Euler = EULER_1D( opts, dimensions, BCs, 'KiU', 1, T_min );

% f_1(u)_x:

f_1 = @( opts, Q, t ) Euler.evaluate_rhs_1( opts, Q, t );

descriptor = {f_1};

term_u_1 = TERM( u_1, {u_1,u_2,u_3}, descriptor, false );

equation_u_1 = EQUATION( u_1, {term_u_1}, 'evolution', '' );

% f_2(u)_x:

f_2 = @( opts, Q, t ) Euler.evaluate_rhs_2( opts, Q, t );

descriptor = {f_2};

term_u_2 = TERM( u_2, {u_1,u_2,u_3}, descriptor, false );

equation_u_2 = EQUATION( u_2, {term_u_2}, 'evolution', '' );

% f_3(u)_x:

f_3 = @( opts, Q, t ) Euler.evaluate_rhs_3( opts, Q, t );

descriptor = {f_3};

term_u_3 = TERM( u_3, {u_1,u_2,u_3}, descriptor, false );

equation_u_3 = EQUATION( u_3, {term_u_3}, 'evolution', '' );

    function [ dt ] = set_dt( pde_system, CFL )
        lev  = pde_system.opts.lev;
        deg  = pde_system.opts.deg;
        dim  = pde_system.unknowns{1}.dimensions{1};
        xMax = dim.max;
        xMin = dim.min;
        dx   = ( xMax - xMin ) / ( 2^lev );
        dt   = CFL * dx / ( 2 * deg - 1 );
    end

%% Create the PDE System:

pde_system = PDE_SYSTEM( opts, {equation_u_1,equation_u_2,equation_u_3}, @set_dt );

end

