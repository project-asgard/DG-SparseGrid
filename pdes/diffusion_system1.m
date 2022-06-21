function pde_system = diffusion_system1( opts )

%
%% Solving the heat equation u_t - u_xx = 0 as a system:
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
soln_t = @(t,p)   exp(-pi^2*t);
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

%% Second unknown (q):

%% Define the analytic solution (optional).

soln_x = @(x,p,t) cos(pi*x)*pi;
soln_t = @(t,p)   exp(-pi^2*t);
soln   = new_md_func(num_dims,{soln_x,soln_t});

%%  Define the boundary conditions

BCL = {soln};
BCR = {soln};

%% Initial conditions

initial_conditions = {soln};

%% Analytic solution

analytic_solution = {soln};

%% Solution Vector

q = UNKNOWN( opts, dimensions, analytic_solution, initial_conditions );

%% Define the terms of the PDE system

dV = @(x,p,t,dat) x.*0+1;
gp = @(x,p,t,dat) x.*0+1;

% d_x q:

pterm_x   = DIV(num_dims,gp,'',0,'N','N','','','',dV);
sd_term_x = SD_TERM({pterm_x});
md_term_x = MD_TERM(num_dims,{sd_term_x});

descriptor = {md_term_x};

term_u = TERM( u, {q}, descriptor, true, false );

equation_u = EQUATION( u, {term_u}, 'evolution', '' );

% d_x u

pterm_x   = GRAD(num_dims,gp,'',0,'D','D','','','',dV);
sd_term_x = SD_TERM({pterm_x});
md_term_x = MD_TERM(num_dims,{sd_term_x});

descriptor = {md_term_x};

term_q = TERM( q, {u}, descriptor, true, false );

equation_q = EQUATION( q, {term_q}, 'closure', '' );

    function [ dt ] = set_dt( pde_system, CFL )
        lev  = pde_system.opts.lev;
        deg  = pde_system.opts.deg;
        dim  = pde_system.unknowns{1}.dimensions{1};
        xMax = dim.max;
        xMin = dim.min;
        dx   = ( xMax - xMin ) / ( 2^lev );
        dt   = CFL * dx^2 / ( 2 * deg - 1 );
    end

pde_system = PDE_SYSTEM( opts, {equation_u,equation_q}, @set_dt );

%%% Hack to plot initial condition %%%
% 
% t=0.0;
% 
% for d=1:num_dims
%     num_fixed_grid = 257;
%     nodes_nodups{d}...
%         = linspace( pde_system.unknowns{1}.dimensions{d}.min,...
%                     pde_system.unknowns{1}.dimensions{d}.max,...
%                     num_fixed_grid );
%     [ Meval{d}, nodes{d}, nodes_count{d} ]...
%         = matrix_plot_D( pde_system.unknowns{1}, pde_system.opts,...
%                          pde_system.unknowns{1}.dimensions{d},...
%                          nodes_nodups{d} );
% end
% 
% x_A  = nodes{1}';
% u0_A = u.analytic_solutions{1}{1}(x_A,1.0,0.0).*u.analytic_solutions{1}{2}(0.0,1);
% uf_A = u.analytic_solutions{1}{1}(x_A,1.0,  t).*u.analytic_solutions{1}{2}(  t,1);
% q0_A = q.analytic_solutions{1}{1}(x_A,1.0,0.0).*q.analytic_solutions{1}{2}(0.0,1);
% qf_A = q.analytic_solutions{1}{1}(x_A,1.0,  t).*q.analytic_solutions{1}{2}(  t,1);
% 
% lo = pde_system.solution_vector.lbounds(1);
% hi = pde_system.solution_vector.ubounds(1);
% u_rs...
%   = wavelet_to_realspace( pde_system.unknowns{1}, pde_system.opts, Meval,...
%                           pde_system.solution_vector.fvec(lo:hi), pde_system.unknowns{1}.hash_table );
% 
% close all
% 
% fig_1 = figure( 1 );
% 
% subplot(2,1,1)
% 
% plot( x_A     , u0_A, ':k', 'linewidth', 2 ); hold on
% plot( x_A     , uf_A, '-k', 'linewidth', 2 )
% plot( nodes{1}, u_rs, '-r', 'linewidth', 2 )
% title( '$u$', 'interpreter', 'latex' )
% 
% lo = pde_system.solution_vector.lbounds(2);
% hi = pde_system.solution_vector.ubounds(2);
% q_rs...
%   = wavelet_to_realspace( pde_system.unknowns{2}, pde_system.opts, Meval,...
%                           pde_system.solution_vector.fvec(lo:hi), pde_system.unknowns{2}.hash_table );
% 
% subplot(2,1,2)
% 
% plot( x_A     , q0_A, ':k', 'linewidth', 2 ); hold on
% plot( x_A     , qf_A, '-k', 'linewidth', 2 )
% plot( nodes{1}, q_rs, '-r', 'linewidth', 2 )
% title( '$q$', 'interpreter', 'latex' )
% 
% exportgraphics( fig_1, 'diffusion_system1.pdf' )
% 
% norm( u_rs - uf_A )/numel(nodes{1})
% norm( q_rs - qf_A )/numel(nodes{1})

%%% End Hack %%%

end

