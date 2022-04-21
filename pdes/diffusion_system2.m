function pde_system = diffusion_system2()

%
opts = OPTS( {} );
opts.lev=6;
opts.fast_FG_matrix_assembly = true;
%
%% First Component: d_t u = - d_x sigma_x - d_y sigma_y + f(u)

%% Define the number of functions evolved:

num_funcs_1 = 1;

dim_x = DIMENSION(0,1);
dim_x.moment_dV = @(x,p,t,dat) 0*x+1;
dim_y = DIMENSION(0,1);
dim_y.moment_dV = @(x,p,t,dat) 0*x+1;
dimensions_1 = {dim_x,dim_y};
num_dims_1 = numel(dimensions_1);

%% Define the analytic solution (optional).

soln_x_1 = @(x,p,t) cos(pi*x)*pi^2;
soln_y_1 = @(y,p,t) 0*y+1;%cos(pi*y);
soln_t_1 = @(t,p)   exp(-2*pi^2*t);
soln_1_1 = new_md_func(num_dims_1,{soln_x_1,soln_y_1,soln_t_1});

soln_1 = {soln_1_1};

%%  Define the boundary conditions

BCL_1 = {soln_1_1};
BCR_1 = {soln_1_1};

%% Initial conditions

initial_conditions_1 = {soln_1_1};

%% Solution Vector

u = UNKNOWN( opts, dimensions_1, num_funcs_1, soln_1, initial_conditions_1 );
u.set_initial_conditions( opts );

%
%% Second Component: (sigma_x,sigma_y) = - (d_x u,d_y u)

%% Define the number of functions evolved:

num_funcs_2 = 2;

%% Define the dimensionality:

dim_x = DIMENSION(0,1);
dim_x.moment_dV = @(x,p,t,dat) 0*x+1;
dim_y = DIMENSION(0,1);
dim_y.moment_dV = @(x,p,t,dat) 0*x+1;
dimensions_2 = {dim_x,dim_y};
num_dims_2 = numel(dimensions_2);

%% Define the analytic solution (optional).

soln_x_1 = @(x,p,t) sin(pi*x)*pi;
soln_y_1 = @(y,p,t) 0*y+1;%cos(pi*y);
soln_t_1 = @(t,p)   exp(-2*pi^2*t);
soln_2_1 = new_md_func(num_dims_2,{soln_x_1,soln_y_1,soln_t_1});

soln_x_2 = @(x,p,t) cos(pi*x);
soln_y_2 = @(y,p,t) sin(pi*y)*pi;
soln_t_2 = @(t,p)   exp(-2*pi^2*t);
soln_2_2 = new_md_func(num_dims_2,{soln_x_2,soln_y_2,soln_t_2});

soln_2 = {soln_2_1,soln_2_2};

%%  Define the boundary conditions

BCL_2 = {soln_2_1,soln_2_2};
BCR_2 = {soln_2_1,soln_2_2};

%% Initial conditions

initial_conditions_2 = {soln_2_1,soln_2_2};

%% Solution Vector

sigma = UNKNOWN( opts, dimensions_2, num_funcs_2, soln_2, initial_conditions_2 );
sigma.set_initial_conditions( opts );

%% Global Solution Vector

Q = {u,sigma};

%% Define the terms of the PDE system

dV = @(x,p,t,dat) 0*x+1;

gm = @(x,p,t,dat) x.*0-1;
gp = @(x,p,t,dat) x.*0+1;

% d_x sigma_x:

subterm_11 = 'ZERO';

pterm_x = DIV(num_dims_2,gp,'',-1,'N','N','','','',dV); % num_dims_2 because it is applied to sigma (think about this?)
pterm_y = MASS(gp);
term_x  = SD_TERM({pterm_x}); % Maybe remove?
term_y  = SD_TERM({pterm_y}); % Maybe remove?
subterm_21 = MD_TERM(num_dims_2,{term_x,term_y});
subterm_22 = 'ZERO';

descriptor = {{subterm_11},{subterm_21,subterm_22}};

term_1 = TERM( u, Q, descriptor );

out = term_1.driver( opts, 0.0 );% Hack to test driver

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

