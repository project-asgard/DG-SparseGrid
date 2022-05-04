function pde_system = diffusion_system1()

%
opts = OPTS( {} );
opts.lev=8;
opts.grid_type = 'FG';
%
%% 

%% Define the number of functions evolved:

dim_x = DIMENSION(0,1);
dim_x.moment_dV = @(x,p,t,dat) 0*x+1;
dimensions = {dim_x};
num_dims = numel(dimensions);

%% First unknown:

%% Define the analytic solution (optional).

soln_x = @(x,p,t) -sin(pi*x)*pi^2;
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

%% Second unknown:

%% Define the analytic solution (optional).

soln_x = @(x,p,t) cos(pi*x)*pi;
soln_t = @(t,p)   exp(-2*pi^2*t);
soln = new_md_func(num_dims,{soln_x,soln_t});

%%  Define the boundary conditions

BCL = {soln};
BCR = {soln};

%% Initial conditions

initial_conditions = {soln};

%% Solution Vector

q = UNKNOWN( opts, dimensions, soln, initial_conditions );
q.set_initial_conditions( opts );

%% Global Solution Vector

Q = {u,q};

%% Define the terms of the PDE system

dV = @(x,p,t,dat) 0*x+1;
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

% d_y sigma_y:

subterm_11 = 'ZERO';

pterm_x = MASS(gp);
pterm_y = DIV(num_dims_2,gp,'',+1,'N','N','','','',dV); % num_dims_2 because it is applied to sigma (think about this?)
term_x  = SD_TERM({pterm_x}); % Maybe remove?
term_y  = SD_TERM({pterm_y}); % Maybe remove?
subterm_21 = 'ZERO';
subterm_22 = MD_TERM(num_dims_2,{term_x,term_y});

descriptor = {{subterm_11},{subterm_21,subterm_22}};

term_2 = TERM( u, Q, descriptor );

% --- Create Equation 1 with term_1 and term_2 here?

% grad_x u

pterm_x = GRAD(num_dims_1,gp,'',-1,'D','D',BCL_2{1},BCR_2{1},'',dV);
pterm_y = MASS(gp);
term_x  = SD_TERM({pterm_x});
term_y  = SD_TERM({pterm_y});
subterm_11 = MD_TERM(num_dims_2,{term_x,term_y});

% mass sigma_x

pterm_x = MASS(gp);
pterm_y = MASS(gp);
term_x  = SD_TERM({pterm_x});
term_y  = SD_TERM({pterm_y});
subterm_21 = MD_TERM(num_dims_2,{term_x,term_y});
subterm_22 = 'ZERO';

descriptor = {{subterm_11},{subterm_21,subterm_22}};

% grad_y u

pterm_x = MASS(gp);
pterm_y = GRAD(num_dims_1,gp,'',-1,'D','D',BCL_2{2},BCR_2{2},'',dV);
term_x  = SD_TERM({pterm_x});
term_y  = SD_TERM({pterm_y});
subterm_11 = MD_TERM(num_dims_2,{term_x,term_y});

% mass sigma_y

pterm_x = MASS(gp);
pterm_y = MASS(gp);
term_x  = SD_TERM({pterm_x});
term_y  = SD_TERM({pterm_y});
subterm_21 = 'ZERO';
subterm_22 = MD_TERM(num_dims_2,{term_x,term_y});

descriptor = {{subterm_11},{subterm_21,subterm_22}};

% assemble PDE terms data

data = {data_1,data_2,data_3,data_4,data_5,data_6};

pde_terms = PDE_TERMS( data );

pde_system = PDE_SYSTEM( {}, {u,sigma}, pde_terms );

pde_system.initialize()

%%% Hack to plot initial condition %%%

for d=1:num_dims_1
    num_fixed_grid = 51;
    nodes_nodups{d}...
        = linspace( pde_system.pde_solutions{1}.dimensions{d}.min,...
                    pde_system.pde_solutions{1}.dimensions{d}.max,...
                    num_fixed_grid );
    [ Meval{d}, nodes{d}, nodes_count{d} ]...
        = matrix_plot_D( pde_system.pde_solutions{1}, pde_system.opts,...
                         pde_system.pde_solutions{1}.dimensions{d},...
                         nodes_nodups{d} );
end

fig_1 = figure( 1 );

u_rs...
  = wavelet_to_realspace( pde_system.pde_solutions{1}, pde_system.opts, Meval,...
                          pde_system.pde_solutions{1}.fval, pde_system.hash_table );

u_rs_nD = singleD_to_multiD(num_dims_1,u_rs,nodes);

subplot(2,2,1)

imagesc( nodes{1}, nodes{2}, u_rs_nD );axis square
title( '$u_0$', 'interpreter', 'latex' )

sigma_x_rs...
  = wavelet_to_realspace( pde_system.pde_solutions{2}, pde_system.opts, Meval,...
                          pde_system.pde_solutions{2}.fval(:,1), pde_system.hash_table );

sigma_x_rs_nD = singleD_to_multiD(num_dims_1,sigma_x_rs,nodes);

subplot(2,2,2)

imagesc( nodes{1}, nodes{2}, sigma_x_rs_nD );axis square
title( '$\sigma_{x,0}$', 'interpreter', 'latex' )

sigma_y_rs...
  = wavelet_to_realspace( pde_system.pde_solutions{2}, pde_system.opts, Meval,...
                          pde_system.pde_solutions{2}.fval(:,2), pde_system.hash_table );

sigma_y_rs_nD = singleD_to_multiD(num_dims_1,sigma_y_rs,nodes);

subplot(2,2,3)

imagesc( nodes{1}, nodes{2}, sigma_y_rs_nD );axis square
title( '$\sigma_{y,0}$', 'interpreter', 'latex' )

sigma_x_tmp = zeros(52,52);
for j = 1 : 52
for i = 1 : 52
    sigma_x_tmp(i,j) = sin(pi*nodes{1}(i))*cos(pi*nodes{2}(j))*pi;
end
end

subplot(2,2,4)

imagesc( nodes{1}, nodes{2}, sigma_x_tmp' );axis square

%%% End Hack %%%

pde_system.pde_solutions = pde_terms.pde_driver( pde_system.pde_solutions );% Remove Later

end

