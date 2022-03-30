function pde_system = diffusion_system2()

%
opts = OPTS( {} );

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

soln_x_1 = @(x,p,t) cos(pi*x);
soln_y_1 = @(y,p,t) cos(pi*y);
soln_t_1 = @(t,p)   exp(-2*pi^2*t);
soln_1_1 = new_md_func(num_dims_1,{soln_x_1,soln_y_1,soln_t_1});

soln_1 = {soln_1_1};

%%  Define the boundary conditions

BCL_1 = {soln_1_1};
BCR_1 = {soln_1_1};

%% Initial conditions

initial_conditions_1 = {soln_1_1};

%% Solution Vector

u = PDE_SOLUTIONS( opts, dimensions_1, soln_1, initial_conditions_1 );

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
soln_y_1 = @(y,p,t) cos(pi*y);
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

sigma = PDE_SOLUTIONS( opts, dimensions_2, soln_2, initial_conditions_2 );

%% Define the terms of the PDE system

dV = @(x,p,t,dat) 0*x+1;

gm = @(x,p,t,dat) x.*0-1;
gp = @(x,p,t,dat) x.*0+1;

% d_x sigma_x:

pterm_x = DIV(num_dims_2,gp,'',+1,'N','N','','','',dV); % num_dims_2 because it is applied to sigma (think about this?)
pterm_y = MASS(gp);
term_x  = SD_TERM({pterm_x}); % Maybe remove?
term_y  = SD_TERM({pterm_y}); % Maybe remove?
data_1.term      = MD_TERM(num_dims_2,{term_x,term_y});
data_1.input_id  = [2,1]; %Meaning apply term to system 2, component 1 [system,component]
data_1.output_id = [1,1]; %Meaning put result on System 1, component 1 [system,component]

% d_y sigma_y:

pterm_x = MASS(gp);
pterm_y = DIV(num_dims_2,gp,'',+1,'N','N','','','',dV); % num_dims_2 because it is applied to sigma (think about this?)
term_x  = SD_TERM({pterm_x}); % Maybe remove?
term_y  = SD_TERM({pterm_y}); % Maybe remove?
data_2.term     = MD_TERM(num_dims_2,{term_x,term_y});
data_2.input_id  = [2,2]; %Meaning apply term to system 2, component 2
data_2.output_id = [1,1]; %Meaning put result on System 1, component 1

% grad_x u

pterm_x = GRAD(num_dims_1,gp,'',-1,'D','D',BCL_2{1},BCR_2{1},'',dV);
pterm_y = MASS(gp);
term_x  = SD_TERM({pterm_x});
term_y  = SD_TERM({pterm_y});
data_3.term = MD_TERM(num_dims_2,{term_x,term_y});
data_3.input_id = [1,1];
data_3.imput_id = [2,1];

% mass sigma_x

pterm_x = MASS(gp);
pterm_y = MASS(gp);
term_x  = SD_TERM({pterm_x});
term_y  = SD_TERM({pterm_y});
data_4.term = MD_TERM(num_dims_2,{term_x,term_y});
data_4.input_id  = [2,1];
data_4.output_id = [2,1];

% grad_y u

pterm_x = MASS(gp);
pterm_y = GRAD(num_dims_1,gp,'',-1,'D','D',BCL_2{2},BCR_2{2},'',dV);
term_x  = SD_TERM({pterm_x});
term_y  = SD_TERM({pterm_y});
data_5.term = MD_TERM(num_dims_2,{term_x,term_y});
data_5.input_id  = [1,1];
data_5.output_id = [2,2];

% mass sigma_y

pterm_x = MASS(gp);
pterm_y = MASS(gp);
term_x  = SD_TERM({pterm_x});
term_y  = SD_TERM({pterm_y});
data_6.term = MD_TERM(num_dims_2,{term_x,term_y});
data_6.input_id  = [2,2];
data_6.output_id = [2,2];

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

