function [err,fval,fval_realspace,nodes,err_realspace,outputs,opts] = ...
    asgard_run_pde(opts,pde)

root_directory = get_root_folder();

tic
figs = [];

num_dims = numel(pde.dimensions);

%% Reset any persistent variables
if opts.time_independent_A || opts.time_independent_build_A || strcmp(opts.timestep_method,'IMEX')
    clear time_advance
end

%% Reset fast multiply
clear fast_2d_matrix_apply

%% Set time step.
dt = pde.set_dt(pde,opts.CFL);
if opts.dt_set_at_runtime
    dt = opts.dt;
else
    opts.dt = dt;
end
if ~opts.quiet; disp(sprintf('dt = %g', dt )); end


%% Construct the Element (Hash) table.
if ~opts.quiet; disp('Constructing hash and inverse hash tables'); end

if opts.use_oldhash
    [HASH,hash_table] = hash_table_nD(pde.get_lev_vec, opts.grid_type);
    pde.hash_table = hash_table;
else
    [elements, elements_idx]    = hash_table_sparse_nD (pde.get_lev_vec, opts.max_lev, opts.grid_type);
    hash_table.elements         = elements;
    hash_table.elements_idx     = elements_idx; % only to get the same order as the old hash table
end


%% (Do not) Construct the connectivity.
if opts.use_connectivity
    pde.connectivity = connect_nD(num_dims,HASH,hash_table,max(pde.get_lev_vec),max(pde.get_lev_vec),opts.grid_type);
else
    connectivity = [];
end

%% Generate 1D mass matrices in each dimension (used in L2 projections)
if ~opts.quiet; disp('Calculate 1D mass matricies in each dimesion'); end
pde = compute_dimension_mass_mat(opts,pde);

%% Generate moment vectors
if ~opts.quiet; disp('Generate moment vector data'); end
pde = calculate_moment_data(pde,opts);

%% Generate initial conditions (both 1D and multi-D).
if ~opts.quiet; disp('Calculate 2D initial condition on the sparse-grid'); end
t = 0;
fval = initial_condition_vector(pde, opts, hash_table, t);
if opts.save_output; fval_t{1} = fval; end

%% Construct the time-independent coefficient matrices
if ~opts.quiet; disp('Calculate time independent matrix coefficients'); end
TD = 0;
pde = get_coeff_mats(pde,opts,t,TD);

%% Construct A_encode / A_data time independent data structures.
if ~opts.quiet; disp('Generate A_encode data structure for time independent coefficients'); end
A_data = global_matrix(pde,opts,hash_table);

%% Create 1D pde file 

if opts.update_params_each_timestep
    
     
    
end

%% Construct Poisson matrix
if ~opts.quiet; disp('Construct matrix for Poisson solve'); end
if pde.solvePoisson
    if DimX>1
        % Construct DeltaX for DimX
    else
        A_Poisson = DeltaX;
    end
end

%% Construct transforms back to realspace for plotting

if opts.build_realspace_output
    for d=1:num_dims
        if strcmp(opts.output_grid,'fixed')
            if d==1
                num_fixed_grid = 51;
            else
                num_fixed_grid = 21;
            end
            nodes_nodups{d} = ...
                linspace(pde.dimensions{d}.min,pde.dimensions{d}.max,num_fixed_grid);
            [Meval{d},nodes{d},nodes_count{d}] = ...
                matrix_plot_D(pde,opts,pde.dimensions{d},nodes_nodups{d});
        elseif strcmp(opts.output_grid,'elements')
            [element_coordinates,element_coordinates_deg] = get_sparse_grid_coordinates(pde,opts,hash_table);
            nodes_nodups{d} = unique(sort(element_coordinates_deg(d,:)));
            [Meval{d},nodes{d},nodes_count{d}] = ...
                matrix_plot_D(pde,opts,pde.dimensions{d},nodes_nodups{d});
        else
            [Meval{d},nodes{d}] = matrix_plot_D(pde,opts,pde.dimensions{d});
            nodes_nodups{d} = nodes{d};
            nodes_count{d} = nodes{d}.*0+1;
        end
    end
end

%% Construct a n-D coordinate array
if opts.build_realspace_output
    coord = get_realspace_coords(pde,nodes);
    coord_nodups = get_realspace_coords(pde,nodes_nodups);
end

%% Plot initial condition
if num_dims <=3
    
    %%
    % Get the real space solution
    
    if opts.build_realspace_output
        fval_realspace = wavelet_to_realspace(pde,opts,Meval,fval,hash_table);
        if ~isempty(pde.solutions)
            fval_realspace_analytic = get_analytic_realspace_solution_D(pde,opts,coord,t);
            fval_realspace_analytic = reshape(fval_realspace_analytic, length(fval_realspace),1);
        end
    end
    
    % construct the moment function handle list for calculating the mass
    if opts.calculate_mass
        [pde,mass,mass_analytic] = calculate_mass(pde,opts,fval,hash_table,t);
        mass_t(1) = mass(1);
    end
    
    if opts.normalize_by_mass && ~isempty(pde.solutions)
        pde.params.norm_fac = mass / mass_analytic(1);
        fval_realspace_analytic = get_analytic_realspace_solution_D(pde,opts,coord,t);
        fval_realspace_analytic = reshape(fval_realspace_analytic, length(fval_realspace), 1);
    end
    
    if opts.build_realspace_output
        f_realspace_nD = singleD_to_multiD(num_dims,fval_realspace,nodes);
        if strcmp(opts.output_grid,'fixed') || strcmp(opts.output_grid,'elements')
            f_realspace_nD = ...
                remove_duplicates(num_dims,f_realspace_nD,nodes_nodups,nodes_count);
        end

        if isempty(pde.solutions)
            f_realspace_analytic_nD = [];
        else
            f_realspace_analytic_nD = get_analytic_realspace_solution_D(pde,opts,coord_nodups,t);
        end
    
        if opts.save_output
            if num_dims <= 3
                f_realspace_nD_t{1} = f_realspace_nD;
            else
                error('Save output for num_dimensions >3 not yet implemented');
            end
        end

        if opts.use_oldhash
        else
            element_coordinates = get_sparse_grid_coordinates(pde,opts,hash_table);
        end

        if norm(fval_realspace) > 0 && ~opts.quiet
            figs.ic = figure('Name','Initial Condition','Units','normalized','Position',[0.7,0.1,0.3,0.3]);

            if opts.use_oldhash
                plot_fval(pde,nodes_nodups,f_realspace_nD,f_realspace_analytic_nD);
            else
                plot_fval(pde,nodes_nodups,f_realspace_nD,f_realspace_analytic_nD,element_coordinates);
            end
        end
    end
    
    %     fval_realspace_SG = real_space_solution_at_coordinates_irregular(pde,fval,coordinates);
    
    if opts.calculate_mass
        [pde,mass,mass_analytic] = calculate_mass(pde,opts,fval,hash_table,t);
        mass_t(1) = mass(1);
    end
    
end

%% Write the initial condition to file.
write_fval = 0;
if write_fval; write_fval_to_file(fval,lev,deg,0); end

%% Check to see if initial resolution meets requested accuracy
if opts.adapt
    if ~opts.quiet
        figs.adapt = figure();
    end
    if opts.adapt_initial_condition
        if ~opts.quiet; disp('Adapting initial for requested accuracy ...'); end
        
        keep_adapting_initial_condition = true;
        while keep_adapting_initial_condition
            num_pre_adapt = numel(fval);
            % first refine
            [pde,fval_tmp,hash_table,A_data,Meval,nodes,nodes_nodups,nodes_count,coord,coord_nodups,~,fval_realspace] ...
                = adapt(pde,opts,figs,fval,hash_table,Meval,nodes,nodes_nodups,nodes_count, ...
                fval_realspace,0,1);
            if num_pre_adapt == numel(fval_tmp)
                keep_adapting_initial_condition = false;
            end
            clear fval_tmp;
            % reproject initial condition onto refined basis
            fval = initial_condition_vector(pde, opts, hash_table, t);
        end
        
        % coarsen
        [pde,~,hash_table,A_data,Meval,nodes,nodes_nodups,nodes_count,coord,coord_nodups] ...
            = adapt(pde,opts,figs,fval,hash_table,Meval,nodes,nodes_nodups,nodes_count, ...
            fval_realspace,1,0);
        % reproject onto coarsend basis
        fval = initial_condition_vector(pde, opts, hash_table, t);
        
    else
        if ~opts.quiet; disp('Checking if initial grid is sufficient for requested accuracy ...'); end
        
        % Check to ensure refinement is not required to start
        pre_refinement_num_DOF = length(fval);
        [~,fval_check] ...
            = adapt(pde,opts,figs,fval,hash_table,Meval,nodes,nodes_nodups,nodes_count,fval_realspace,0,1);
        if (length(fval_check)>pre_refinement_num_DOF)
            %             error('Initial grid was insifficient for requested accuracy');
        end
        clear fval_check;
    end
else
    if ~opts.quiet; disp(['Number of DOF : ', num2str(numel(fval))]); end
end

fval_analytic = exact_solution_vector(pde,opts,hash_table,t);
err_wavelet = sqrt(mean((fval(:) - fval_analytic(:)).^2));
if ~opts.quiet
    disp(['    num_dof : ', num2str(numel(fval))]);
    disp(['    wavelet space absolute err : ', num2str(err_wavelet)]);
    disp(['    wavelet space relative err : ', num2str(err_wavelet/max(abs(fval_analytic(:)))*100), ' %']);
    disp(['    wavelet space absolute err (2-norm) : ', num2str(norm(fval-fval_analytic))]);
end

if num_dims <=3 && opts.build_realspace_output
    fval_realspace = wavelet_to_realspace(pde,opts,Meval,fval,hash_table);
    fval_realspace_analytic = get_analytic_realspace_solution_D(pde,opts,coord,t);
    err_realspace = sqrt(mean((fval_realspace(:) - fval_realspace_analytic(:)).^2));
    if ~opts.quiet
        disp(['    real space absolute err : ', num2str(err_realspace)]);
        disp(['    real space relative err : ', num2str(err_realspace/max(abs(fval_realspace_analytic(:)))*100), ' %']);
        disp(['    real space absolute err (2-norm) : ', num2str(norm(fval_realspace(:)-fval_realspace_analytic(:)))]);
    end
end

% need to clean up this interface!
if opts.build_realspace_output
    outputs = save_output([],0,pde,opts,num_dims,fval,fval_realspace,f_realspace_analytic_nD,nodes,nodes_nodups,nodes_count,t,dt,toc,root_directory,hash_table);
else
    outputs.fval = fval;
end
if opts.calculate_mass
    outputs.mass_t = mass_t;
end

%% Time Loop
count=1;
err = 1e9;
sol = 0;
if opts.record_video
    v = VideoWriter('video.avi','MPEG-4');
    open(v);
end
if ~opts.quiet; disp('Advancing time ...'); end
for L = 1:opts.num_steps
    
    tic;
    if ~opts.quiet
        timeStr = sprintf('Step %i of %i (t= %f) (dof=%i)',L,opts.num_steps,t,numel(fval));
        disp(timeStr);
    end
    
    Emax = 0;
    
    % Coarsen Grid
    if opts.adapt
        [pde,fval,hash_table,A_data,Meval,nodes,nodes_nodups,nodes_count,coord,coord_nodups] ...
            = adapt(pde,opts,figs,fval,hash_table,Meval,nodes,nodes_nodups,nodes_count,fval_realspace,1,0);
    end
    
    needs_adapting = true;
    while needs_adapting
        %     for adapt_step = 1:2
        
        if pde.solvePoisson
            %%
            % Solve Poisson to get E (from 1-rho=1-int f dv)
            if ~quiet; disp('    Solve poisson to get E'); end
            %[E,u] = PoissonSolve2(LevX,Deg,Lmax,fval,A_Poisson,FMWT_COMP_x,Vmax,index1D);
            [E,u] = PoissonSolve(LevX,deg,Lmax,fval,A_Poisson,FMWT_COMP_x,Vmax);
            Emax = max(abs(Meval{2}*E)); % TODO : this clearly is problem dependent
        end
        
        if pde.applySpecifiedE
            %%
            % Apply specified E
            if ~quiet; disp('    Apply specified E'); end
            E = forwardMWT(LevX,deg,Lmin,Lmax,pde.Ex,pde.params);
            E = E * pde.Et(t,params);
            Emax = max(abs(Meval{2}*E)); % TODO : this clearly is problem dependent
        end
               
        if ~opts.quiet; disp('    Calculate time dependent matrix coeffs'); end
        if num_dims==2
            if (pde.applySpecifiedE || pde.solvePoisson)
                
                %%
                % Generate EMassX time dependent coefficient matrix.
                
                EMassX = matrix_coeff_TD(LevX,deg,Lmin,Lmax,E,pde.transform_blocks);
                
                %%
                % Set the dat portion of the EMassX part of E.d_dv term.
                
                pde.terms{2}{2}.dat = E;
                
            end
        end
        
        %%
        % Now construct the TD coeff_mats.
        
        TD = 1;
        pde = get_coeff_mats(pde,opts,t,TD);
        
        %%
        % Advance in time
        
        if ~opts.quiet; disp(['    Time step (dt=',num2str(dt),')']); end
        
        %%
        % Write the A_data structure components for use in HPC version.
        write_A_data = 0;
        if write_A_data && L==1; write_A_data_to_file(A_data,lev,deg); end
        
        fval_unstepped = fval;
        if strcmp(opts.timestep_method,'ode15s') && L>1
            fval = deval(sol,t+dt);
        else
            [fval,sol] = time_advance(pde,opts,A_data,fval,t,dt,opts.deg,hash_table,[],[]);
        end

        if opts.build_realspace_output
            fval_realspace = wavelet_to_realspace(pde,opts,Meval,fval,hash_table);
        end
        
        %%
        % Refine Grid - determine which elements to add, but reset fval to
        % fval_previous with those new elements added and time_advance
        % again
        
        if opts.adapt
            if ~opts.quiet; disp('Adapt grid ...'); end
            
            num_elements_0 = numel(fval);
            
            [pde,~,hash_table,A_data,Meval,nodes,nodes_nodups,nodes_count,coord,coord_nodups,fval_unstepped_adapted] ...
                = adapt(pde,opts,figs,fval,hash_table,Meval,nodes,nodes_nodups,nodes_count, ...
                fval_realspace,0,1,fval_unstepped);
            
            for d=1:num_dims
                assert(numel(nodes_nodups{d})==numel(nodes_count{d}))
            end
            
            num_elements_adapted = numel(fval_unstepped_adapted);
            
            if num_elements_0 == num_elements_adapted
                if ~opts.quiet
                    disp(['No more adaption needed - advancing time']);
                    disp(['    t = ', num2str(t)]);
                end
                needs_adapting = false;
                fval = time_advance(pde,opts,A_data,fval_unstepped_adapted,t,dt,opts.deg,hash_table,[],[]);
                if opts.build_realspace_output
                    fval_realspace = wavelet_to_realspace(pde,opts,Meval,fval,hash_table);
                end
            else
                if ~opts.quiet
                    disp(['Still needs adaption iterating ... added ', ...
                        num2str(num_elements_adapted-num_elements_0),' elements'])
                    disp(['    t = ', num2str(t)]);
                end
                fval = fval_unstepped_adapted;
            end
        else
            needs_adapting = false;
        end
        
    end
    
    %%
    % Write the present fval to file.
    if write_fval; write_fval_to_file(fval,lev,deg,L); end
    
    if num_dims <=3 && opts.build_realspace_output
        
        %%
        % Get the real space solution
        
        fval_realspace = wavelet_to_realspace(pde,opts,Meval,fval,hash_table);
        if opts.calculate_mass
            [pde,mass,~] = calculate_mass(pde,opts,fval,hash_table,t);
            mass_t(L+1) = mass(1);
            outputs.mass_t = mass_t;
        end
        
        %%
        % Try with function convertToRealSpace
        
        tryConvertToRealSpace = 0;
        if tryConvertToRealSpace
            LminB = zeros(1,num_dims);
            LmaxB = zeros(1,num_dims);
            for d=1:num_dims
                LminB(d) = pde.dimensions{d}.min;
                LmaxB(d) = pde.dimensions{d}.max;
            end
            fval_realspaceB = convert_to_real_space(pde,num_dims,lev,deg,gridType,LminB,LmaxB,fval,lev);
            % fval_realspace = fval_realspaceB;
        end
        
        if opts.use_oldhash
        else
            element_coordinates = get_sparse_grid_coordinates(pde,opts,hash_table);
        end
        % fval_realspace_SG = real_space_solution_at_coordinates(pde,fval,coordinates);
        
    end
    
    %%
    % Check against known solution
    
    if opts.calculate_mass && ~opts.quiet
        disp(['    total integrated mass : ', num2str(mass(1))]);
    end
    if ~isempty(pde.solutions)
        
        %%
        % Check the wavelet space solution with the analytic solution
        
        fval_analytic = exact_solution_vector(pde,opts,hash_table,t+dt);
        
        assert(numel(fval)==numel(fval_analytic));
        
        err_wavelet = sqrt(mean((fval(:) - fval_analytic(:)).^2));
        outputs.err{L+1} = err_wavelet;
        outputs.rel_err{L+1} = err_wavelet/norm(fval_analytic(:),Inf);
        if ~opts.quiet  
            disp(['    num_dof : ', num2str(numel(fval))]);
            disp(['    wavelet space absolute err : ', num2str(err_wavelet)]);
            disp(['    wavelet space relative err : ', num2str(err_wavelet/norm(fval_analytic(:),Inf)*100), ' %']);
        end
        
        %%
        % Check the realspace solution
        
        if num_dims <= 3 && opts.build_realspace_output
            if ~opts.quiet
                disp(['t: ',num2str(t)]);
                disp(['dt: ',num2str(dt)]);
            end
            
            fval_realspace_analytic = get_analytic_realspace_solution_D(pde,opts,coord,t+dt);
            err_realspace = sqrt(mean((fval_realspace(:) - fval_realspace_analytic(:)).^2));
            if ~opts.quiet
                disp(['    real space absolute err : ', num2str(err_realspace)]);
                disp(['    real space relative err : ', num2str(err_realspace/max(abs(fval_realspace_analytic(:)))*100), ' %']);
            end
        end
        
        catch_min_error = false;
        if catch_min_error && err < err_wavelet
            disp('Error is now going up?');
        end
        err = err_wavelet;
    end
    
    % Reshape realspace solution and plot
    
    if num_dims <= 3 && opts.build_realspace_output
        
        f_realspace_nD = singleD_to_multiD(num_dims,fval_realspace,nodes);
        if strcmp(opts.output_grid,'fixed') || strcmp(opts.output_grid,'elements')
            f_realspace_nD = ...
                remove_duplicates(num_dims,f_realspace_nD,nodes_nodups,nodes_count);
        end
        
        f_realspace_analytic_nD = get_analytic_realspace_solution_D(pde,opts,coord_nodups,t+dt);
        
        element_coordinates = [];
        if opts.use_oldhash
        else
            element_coordinates = get_sparse_grid_coordinates(pde,opts,hash_table);
        end
        
        %%
        % Plot results
        
        if mod(L,opts.plot_freq)==0 && ~opts.quiet
            
            if isfield(figs,'solution')
                figure(figs.solution);
            else
                %figs.solution = figure('Name','Solution','Units','normalized','Position',[0.1,0.1,0.5,0.5]);
                figs.solution = figure('Name','Solution','Units','normalized','Position',[0.1,0.1,0.3,0.3]);
            end
                       
            plot_fval(pde,nodes_nodups,f_realspace_nD,f_realspace_analytic_nD,element_coordinates);
            
            if opts.record_video
                frame = getframe(gcf);
                writeVideo(v,frame);
            end
            
            % this is just for the RE paper
            plot_fval_in_cyl = false;
            if plot_fval_in_cyl
                p = nodes{1};
                z = nodes{2};
                f = reshape(fval_realspace,size(fval_realspace_analytic));
                pper = linspace(0,max(p),100);
                ppar = linspace(-max(p),+max(p),201);
                [ppar2d,pper2d] = meshgrid(ppar,pper);
                p2dA = sqrt(ppar2d.^2+pper2d.^2);
                z2dA = cos(atan2(pper2d,ppar2d));
                f2d = interp2(p,z,f,p2dA,z2dA,'spline',0);
                levs = linspace(1,10,10)./10.*max(f(:));
                figure(87)
                contour(ppar,pper,f2d,levs)
            end
                   
        end       
    end
      
    % Hack right now to update the params data every time step
    
    if opts.update_params_each_timestep
        
        %f_p0 = f_realspace_nD(:,1); % get the f(0,z) value (or as close to p=0 as the nodes allow)
        %       alpha_z = @(z) (2/sqrt(pi)-interp1(nodes{2},f_p0,z,'spline','extrap'))/dt;
        %alpha_z = @(z) (pde.params.f0_p(nodes_nodups{1}(1))-interp1(nodes_nodups{2},f_p0,z,'spline','extrap'))/dt;
        %pde.params.alpha_z = alpha_z;
        %z = linspace(-1,1,100);
        %outputs.alpha_t0{L+1} = alpha_z;
        %         outputs.alpha_t{L+1} = sum(alpha_z(z));
        %outputs.alpha_t{L+1} = alpha_z(0);
        
        %Get 1D hash
        hash_table_1D = hash_table_2D_to_1D(hash_table,opts);
        
        %Make fake pde file
        pde_1d.dimensions = pde.dimensions(1);
        
        
        mom0 = moment_mat{1}*fval; %integral of (f,1)_v
        mom0_real = wavelet_to_realspace(pde_1d,opts,Meval(1),mom0,hash_table_1D);
        %mom0_vals = singleD_to_multiD(num_dims-1,mom0_real,nodes(1));
        pde.params.n  = @(x) interp1(nodes{1},mom0_real,x,'nearest','extrap');
        
        mom1 = moment_mat{2}*fval; %integral of (f,v)_v
        mom1_real = wavelet_to_realspace(pde_1d,opts,Meval(1),mom1,hash_table_1D);
        pde.params.u  = @(x) interp1(nodes{1},mom1_real,x,'nearest','extrap')./pde.params.n(x);
        
        mom2 = moment_mat{3}*fval; %integral of (f,v^2)_v
        mom2_real = wavelet_to_realspace(pde_1d,opts,Meval(1),mom2,hash_table_1D);
        pde.params.th = @(x) interp1(nodes{1},mom2_real,x,'nearest','extrap')./pde.params.n(x) - pde.params.u(x).^2;
        
        figure(1000);
        subplot(2,2,1);
        plot(nodes{1},pde.params.n(nodes{1}));
        title('n_f');
        subplot(2,2,2);
        plot(nodes{1},pde.params.u(nodes{1}));
        title('u_f');
        subplot(2,2,3);
        plot(nodes{1},pde.params.th(nodes{1}));
        title('th_f');
        sgtitle('Fluid Variables');
        
        figure(1001);
        subplot(2,2,1);
        plot(nodes{1},mom0_real);
        title('(f,1)_v');
        subplot(2,2,2);
        plot(nodes{1},mom1_real);
        title('(f,v)_v');
        subplot(2,2,3);
        plot(nodes{1},mom2_real);
        title('(f,v^2)_v');
        sgtitle('Moments');
        
        
        
    end
    
    % Save output
        
    count=count+1;
    t1 = toc;
    if ~opts.quiet; disp(['Took ' num2str(t1) ' [s]']); end
       
    if opts.build_realspace_output
        outputs = save_output(outputs,L,pde,opts,num_dims,fval,fval_realspace,f_realspace_analytic_nD,nodes,nodes_nodups,nodes_count,t,dt,t1,root_directory,hash_table);
    else
        outputs.fval = fval;
    end
    
    t = t + dt;
       
    if opts.calculate_mass && ~opts.quiet
        if isfield(figs,'mass')
            figure(figs.mass);
        else
            figs.mass = figure('Name','Moments','Units','normalized','Position',[0.7,0.5,0.2,0.3]);
        end
        plot_val = get_moments_for_plot(pde);
        semilogy(abs(plot_val));
    end
         
end

if opts.record_video
   close(v); 
end

outputs.pde = pde;
outputs.opts = opts;
if ~opts.build_realspace_output %Set realspace values to empty
    fval_realspace = [];
    nodes = [];
    err_realspace = [];
end

end
