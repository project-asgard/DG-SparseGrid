function pde_system = euler_system1()

%
opts = OPTS( {} );
opts.lev=8;
opts.deg=3;
opts.grid_type = 'FG';
opts.fast_FG_matrix_assembly = true;
opts.timestep_method = 'SSPRK3';
%

TEST = 'Riemann';

assert(any(strcmp(TEST,{'LinearAdvection','Riemann'})))

%
%% Solving the Euler equations as a system:
%  u_t + f(u)_x = 0

%% Define the dimensionality:

dim_x = DIMENSION(0,1);
dim_x.moment_dV = @(x,p,t,dat) 0*x+1;
dimensions = {dim_x};
num_dims = numel(dimensions);

%% First unknown (density: u_1):

%% Define the analytic solution (optional).

switch TEST
    case 'LinearAdvection'
        soln_x = @(x,p,t) 1.+.5*sin(2*pi*x);
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
        soln_x = @(x,p,t) 1.+.5*sin(2*pi*x);
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
        soln_x = @(x,p,t) 0.5*(1.+.5*sin(2*pi.*x)).*(1.+1.e-10);
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

Euler = EULER_1D( opts, dimensions, BCs, 'LLF', 1, T_min );

% f_1(u)_x:

f_1 = @( opts, Q, t, output_unknown ) Euler.evaluate_rhs_1( opts, Q, t, output_unknown );

descriptor = {f_1};

term_u_1 = TERM( u_1, {u_1,u_2,u_3}, descriptor, false );

equation_u_1 = EQUATION( u_1, {term_u_1}, 'evolution', '' );

% f_2(u)_x:

f_2 = @( opts, Q, t, output_unknown ) Euler.evaluate_rhs_2( opts, Q, t, output_unknown );

descriptor = {f_2};

term_u_2 = TERM( u_2, {u_1,u_2,u_3}, descriptor, false );

equation_u_2 = EQUATION( u_2, {term_u_2}, 'evolution', '' );

% f_3(u)_x:

f_3 = @( opts, Q, t, output_unknown ) Euler.evaluate_rhs_3( opts, Q, t, output_unknown );

descriptor = {f_3};

term_u_3 = TERM( u_3, {u_1,u_2,u_3}, descriptor, false );

equation_u_3 = EQUATION( u_3, {term_u_3}, 'evolution', '' );

%% Create the PDE System:

pde_system = PDE_SYSTEM( opts, {equation_u_1,equation_u_2,equation_u_3} );

pde_system.set_initial_conditions;

switch TEST
    case 'LinearAdvection'
        t_f = 1.0;
    case 'Riemann'
        t_f = 0.1;
end

t      = 0.0;
iCycle = 0;
while t < t_f
    
    iCycle = iCycle + 1;
    
    dt = 0.1/((2*opts.deg-1)*(2^opts.lev));
    
    if( t + dt > t_f )
        
        dt = t_f - t;
        
    end
    
    time_stepper( pde_system, t, dt );
    
    t = t + dt;
    
    if mod(iCycle,100)==0; fprintf('t, dt = %f %f\n', t, dt); end
    
end

%%% Hack to plot initial condition %%%

for d=1:num_dims
    num_fixed_grid = 257;
    nodes_nodups{d}...
        = linspace( pde_system.unknowns{1}.dimensions{d}.min,...
                    pde_system.unknowns{1}.dimensions{d}.max,...
                    num_fixed_grid );
    [ Meval{d}, nodes{d}, nodes_count{d} ]...
        = matrix_plot_D( pde_system.unknowns{1}, pde_system.opts,...
                         pde_system.unknowns{1}.dimensions{d},...
                         nodes_nodups{d} );
end

close all

x_A  = nodes{1}';
u_1_A = u_1.analytic_solutions{1}{1}(x_A,1.0,0.0).*u_1.analytic_solutions{1}{2}(0.0,1);
u_2_A = u_2.analytic_solutions{1}{1}(x_A,1.0,0.0).*u_2.analytic_solutions{1}{2}(0.0,1);
u_3_A = u_3.analytic_solutions{1}{1}(x_A,1.0,0.0).*u_3.analytic_solutions{1}{2}(0.0,1);

% --- Exact Solution to Riemann Problem ---

[ x_E, D_E, U_E, T_E, ~ ] = ReadAnalyticRiemann( 'RiemannProblemExact.out' );
u_1_E = D_E;
u_2_E = D_E .* U_E;
u_3_E = 0.5 .* D_E .* ( U_E.^2 + T_E );

lo = pde_system.solution_vector.lbounds(1);
hi = pde_system.solution_vector.ubounds(1);
u_1_N...
  = wavelet_to_realspace( pde_system.unknowns{1}, pde_system.opts, Meval,...
                          pde_system.solution_vector.fvec(lo:hi), pde_system.unknowns{1}.hash_table );
lo = pde_system.solution_vector.lbounds(2);
hi = pde_system.solution_vector.ubounds(2);
u_2_N...
  = wavelet_to_realspace( pde_system.unknowns{2}, pde_system.opts, Meval,...
                          pde_system.solution_vector.fvec(lo:hi), pde_system.unknowns{2}.hash_table );
lo = pde_system.solution_vector.lbounds(3);
hi = pde_system.solution_vector.ubounds(3);
u_3_N...
  = wavelet_to_realspace( pde_system.unknowns{3}, pde_system.opts, Meval,...
                          pde_system.solution_vector.fvec(lo:hi), pde_system.unknowns{3}.hash_table );

fig_1 = figure( 1 );

subplot(3,1,1)

plot( x_A     , u_1_A, ':k', 'linewidth', 2 ); hold on
if(strcmp(TEST,'Riemann'))
plot( x_E     , u_1_E, '-r', 'linewidth', 2 )
end
plot( nodes{1}, u_1_N, '-k', 'linewidth', 2 )
title( '$u_{1}$', 'interpreter', 'latex' )

subplot(3,1,2)

plot( x_A     , u_2_A, ':k', 'linewidth', 2 ); hold on
if(strcmp(TEST,'Riemann'))
plot( x_E     , u_2_E, '-r', 'linewidth', 2 )
end
plot( nodes{1}, u_2_N, '-k', 'linewidth', 2 )
title( '$u_2$', 'interpreter', 'latex' )

subplot(3,1,3)

plot( x_A     , u_3_A, ':k', 'linewidth', 2 ); hold on
if(strcmp(TEST,'Riemann'))
plot( x_E     , u_3_E, '-r', 'linewidth', 2 )
end
plot( nodes{1}, u_3_N, '-k', 'linewidth', 2 )
title( '$u_3$', 'interpreter', 'latex' )

exportgraphics( fig_1, 'euler_system1.pdf' )

format long

[norm( u_1_N - u_1_A )/numel(nodes{1}),...
 norm( u_2_N - u_2_A )/numel(nodes{1}),...
 norm( u_3_N - u_3_A )/numel(nodes{1}),...
 t]

%%% End Hack %%%

end

