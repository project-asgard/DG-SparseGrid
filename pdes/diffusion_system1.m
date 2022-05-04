function pde_system = diffusion_system1()

%
opts = OPTS( {} );
opts.lev=4;
opts.grid_type = 'FG';
%

%
%% Solving the heat equation u_t + u_xx = 0 as a system:
%  u_t + q_x = 0
%    q + u_x = 0

%% Define the dimensionality:

dim_x = DIMENSION(0,1);
dim_x.moment_dV = @(x,p,t,dat) 0*x+1;
dimensions = {dim_x};
num_dims = numel(dimensions);

%% First unknown (u):

%% Define the analytic solution (optional).

soln_x = @(x,p,t) sin(pi*x);
soln_t = @(t,p)   exp(-2*pi^2*t);
soln   = new_md_func(num_dims,{soln_x,soln_t});

%%  Define the boundary conditions

BCL = {soln};
BCR = {soln};

%% Initial conditions

initial_conditions = {soln};

%% Analytic solution

analytic_solution = {soln};

%% Solution Vector

u = UNKNOWN( opts, dimensions, analytic_solution, initial_conditions );
u.set_initial_conditions( opts );

%% Second unknown (q):

%% Define the analytic solution (optional).

soln_x = @(x,p,t) cos(pi*x)*pi;
soln_t = @(t,p)   exp(-2*pi^2*t);
soln = new_md_func(num_dims,{soln_x,soln_t});

%%  Define the boundary conditions

BCL = {soln};
BCR = {soln};

%% Initial conditions

initial_conditions = {soln};

%% Analytic solution

analytic_solution = {soln};

%% Solution Vector

q = UNKNOWN( opts, dimensions, soln, initial_conditions );
q.set_initial_conditions( opts );

%% Global Solution Vector

Q = {u,q};

%% Define the terms of the PDE system

dV = @(x,p,t,dat) x.*0+1;
gp = @(x,p,t,dat) x.*0+1;

% d_x q:

pterm_x   = DIV(num_dims,gp,'',0,'N','N','','','',dV);
sd_term_x = SD_TERM({pterm_x});
md_term_x = MD_TERM(num_dims,{sd_term_x});

descriptor = {md_term_x};

term_u = TERM( u, {q}, descriptor );

equation_u = EQUATION( u, {term_u}, 'evolution', '' );

% d_x u

pterm_x   = GRAD(num_dims,gp,'',0,'N','N','','','',dV);
sd_term_x = SD_TERM({pterm_x});
md_term_x = MD_TERM(num_dims,{sd_term_x});

descriptor = {md_term_x};

term_q = TERM( q, {u}, descriptor );

equation_q = EQUATION( q, {term_q}, 'closure', '' );

pde_system = PDE_SYSTEM( opts, {equation_u,equation_q} );

tic
out = term_u.driver( opts, 0.0 );% Hack to test driver
toc

%%%%  Change opts.lev from 3 to 4 to see drop in error
norm(out-u.fval)

tic
out = equation_u.LHS_term.driver( opts, 0 );
toc

norm(out-u.fval)

%%% Hack to plot initial condition %%%

for d=1:num_dims
    num_fixed_grid = 51;
    nodes_nodups{d}...
        = linspace( pde_system.unknowns{1}.dimensions{d}.min,...
                    pde_system.unknowns{1}.dimensions{d}.max,...
                    num_fixed_grid );
    [ Meval{d}, nodes{d}, nodes_count{d} ]...
        = matrix_plot_D( pde_system.unknowns{1}, pde_system.opts,...
                         pde_system.unknowns{1}.dimensions{d},...
                         nodes_nodups{d} );
end

u_rs...
  = wavelet_to_realspace( pde_system.unknowns{1}, pde_system.opts, Meval,...
                          pde_system.unknowns{1}.fval, pde_system.unknowns{1}.hash_table );

subplot(2,1,1)

plot( nodes{1}, u_rs, '-k', 'linewidth', 2 )
title( '$u_0$', 'interpreter', 'latex' )

q_rs...
  = wavelet_to_realspace( pde_system.unknowns{2}, pde_system.opts, Meval,...
                          pde_system.unknowns{2}.fval, pde_system.unknowns{2}.hash_table );

subplot(2,1,2)

plot( nodes{1}, q_rs, '-k', 'linewidth', 2 )
title( '$q_0$', 'interpreter', 'latex' )

%%% End Hack %%%

end

