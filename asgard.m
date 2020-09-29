%% MATLAB (reference) version of the ASGarD solver

function [err,fval,fval_realspace,nodes,err_realspace] = asgard(pde,varargin)

tic

format short e
folder = fileparts(which(mfilename));
addpath(genpath(folder));
root_directory = get_root_folder();

%% Load PDE and runtime defaults
runtime_defaults

%% Reset any persistent variables
if opts.time_independent_A | opts.time_independent_build_A
    clear time_advance
end

%% Check PDE
pde = check_pde(pde,opts);

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
    [HASH,hash_table] = hash_table_nD(opts.lev_vec, opts.grid_type);
    pde.hash_table = hash_table;
else
    [elements, elements_idx]    = hash_table_sparse_nD (opts.lev_vec, opts.max_lev, opts.grid_type);
    hash_table.elements         = elements;
    hash_table.elements_idx     = elements_idx; % only to get the same order as the old hash table
end

%% Construct the 1D multi-wavelet transform (stored as dense matrix blocks)
[~, pde.transform_blocks] = OperatorTwoScale_wavelet2(opts.deg, opts.max_lev);


%% (Do not) Construct the connectivity.
if opts.use_connectivity
    pde.connectivity = connect_nD(num_dimensions,HASH,hash_table,max(opts.lev_vec),max(opts.lev_vec),opts.grid_type);
else
    connectivity = [];
end


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
for d=1:num_dimensions
    if strcmp(opts.output_grid,'fixed')
        num_fixed_grid = 51;
        nodes_nodups{d} = ...
            linspace(pde.dimensions{d}.domainMin,pde.dimensions{d}.domainMax,num_fixed_grid);
        [Meval{d},nodes{d},nodes_count{d}] = ...
            matrix_plot_D(pde,opts,pde.dimensions{d},nodes_nodups{d});
    else
        [Meval{d},nodes{d}] = matrix_plot_D(pde,opts,pde.dimensions{d});
        nodes_nodups{d} = nodes{d};
        nodes_count{d} = nodes{d}.*0+1;
    end
end

%% Construct a n-D coordinate array
coord = get_realspace_coords(pde,nodes);
coord_nodups = get_realspace_coords(pde,nodes_nodups);


%% Plot initial condition
if num_dimensions <=3
    
    %%
    % Get the real space solution
    fval_realspace = wavelet_to_realspace(pde,opts,Meval,fval,hash_table);    
    f_realspace_nD = singleD_to_multiD(num_dimensions,fval_realspace,nodes);   
    if strcmp(opts.output_grid,'fixed')
        f_realspace_nD = ...
            remove_duplicates(num_dimensions,f_realspace_nD,nodes_nodups,nodes_count);
    end
    
    f_realspace_analytic_nD = get_analytic_realspace_solution_D(pde,opts,coord_nodups,t);

    if opts.save_output
        if num_dimensions <= 3
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
        if opts.use_oldhash
            plot_fval(pde,nodes_nodups,f_realspace_nD,f_realspace_analytic_nD);
        else
            plot_fval(pde,nodes_nodups,f_realspace_nD,f_realspace_analytic_nD,element_coordinates);
        end
    end
    
    %     fval_realspace_SG = real_space_solution_at_coordinates_irregular(pde,fval,coordinates);
    
end

%% Write the initial condition to file.
write_fval = 0;
if write_fval; write_fval_to_file(fval,lev,deg,0); end

%% Check to see if initial resolution meets requested accuracy
if opts.adapt 
    if opts.adapt_initial_condition
        if ~opts.quiet; disp('Adapting initial for requested accuracy ...'); end
                       
        keep_adapting_initial_condition = true;
        while keep_adapting_initial_condition
            num_pre_adapt = numel(fval);
            % first refine
            [pde,fval_tmp,hash_table,A_data,Meval,nodes,nodes_nodups,nodes_count,coord,coord_nodups,~,fval_realspace] ...
                = adapt(pde,opts,fval,hash_table,Meval,nodes,nodes_nodups,nodes_count, ...
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
            = adapt(pde,opts,fval,hash_table,Meval,nodes,nodes_nodups,nodes_count, ...
            fval_realspace,1,0);
        % reproject onto coarsend basis
        fval = initial_condition_vector(pde, opts, hash_table, t);

    else
        if ~opts.quiet; disp('Checking if initial grid is sufficient for requested accuracy ...'); end
        
        % Check to ensure refinement is not required to start
        pre_refinement_num_DOF = length(fval);
        [~,fval_check] ...
            = adapt(pde,opts,fval,hash_table,Meval,nodes,nodes_nodups,nodes_count,fval_realspace,0,1);
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

if num_dimensions <=3
    fval_realspace = wavelet_to_realspace(pde,opts,Meval,fval,hash_table);
    fval_realspace_analytic = get_analytic_realspace_solution_D(pde,opts,coord,t);
    err_realspace = sqrt(mean((fval_realspace(:) - fval_realspace_analytic(:)).^2));
    if ~opts.quiet     
        disp(['    real space absolute err : ', num2str(err_realspace)]);
        disp(['    real space relative err : ', num2str(err_realspace/max(abs(fval_realspace_analytic(:)))*100), ' %']);
        disp(['    real space absolute err (2-norm) : ', num2str(norm(fval_realspace(:)-fval_realspace_analytic(:)))]);
    end
end

%% Time Loop
count=1;
err = 1e9;
if ~opts.quiet; disp('Advancing time ...'); end
for L = 1:num_steps
    
    tic;
    timeStr = sprintf('Step %i of %i at %f seconds',L,num_steps,t);
    
    if ~opts.quiet; disp(timeStr); end
    Emax = 0;
    
    % Coarsen Grid
    if opts.adapt
        [pde,fval,hash_table,A_data,Meval,nodes,nodes_nodups,nodes_count,coord,coord_nodups] ...
            = adapt(pde,opts,fval,hash_table,Meval,nodes,nodes_nodups,nodes_count,fval_realspace,1,0);
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
        if num_dimensions==2
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
        fval = time_advance(pde,opts,A_data,fval,t,dt,opts.deg,hash_table,[],[]);
        fval_realspace = wavelet_to_realspace(pde,opts,Meval,fval,hash_table);
        
        %%
        % Refine Grid - determine which elements to add, but reset fval to
        % fval_previous with those new elements added and time_advance
        % again
        
        if opts.adapt
            if ~opts.quiet; disp('Adapt grid ...'); end
            
            num_elements_0 = numel(fval);
            
            [pde,~,hash_table,A_data,Meval,nodes,nodes_nodups,nodes_count,coord,coord_nodups,fval_unstepped_adapted] ...
                = adapt(pde,opts,fval,hash_table,Meval,nodes,nodes_nodups,nodes_count, ...
                fval_realspace,0,1,fval_unstepped);
            
            num_elements_adapted = numel(fval_unstepped_adapted);
            
            if num_elements_0 == num_elements_adapted
                if ~opts.quiet
                    disp(['No more adaption needed - advancing time']);
                    disp(['    t = ', num2str(t)]);
                end
                needs_adapting = false;
                fval = time_advance(pde,opts,A_data,fval_unstepped_adapted,t,dt,opts.deg,hash_table,[],[]);
                fval_realspace = wavelet_to_realspace(pde,opts,Meval,fval,hash_table);
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
    
    if num_dimensions <=3
        
        %%
        % Get the real space solution
        
        fval_realspace = wavelet_to_realspace(pde,opts,Meval,fval,hash_table);

        %%
        % Try with function convertToRealSpace
        
        tryConvertToRealSpace = 0;
        if tryConvertToRealSpace
            LminB = zeros(1,num_dimensions);
            LmaxB = zeros(1,num_dimensions);
            for d=1:num_dimensions
                LminB(d) = pde.dimensions{d}.domainMin;
                LmaxB(d) = pde.dimensions{d}.domainMax;
            end
            fval_realspaceB = convert_to_real_space(pde,num_dimensions,lev,deg,gridType,LminB,LmaxB,fval,lev);
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
    if pde.checkAnalytic
        
        %%
        % Check the wavelet space solution with the analytic solution
        
        fval_analytic = exact_solution_vector(pde,opts,hash_table,t+dt);
        
        assert(numel(fval)==numel(fval_analytic));
        
        err_wavelet = sqrt(mean((fval(:) - fval_analytic(:)).^2));
        if ~opts.quiet  
            disp(['    num_dof : ', num2str(numel(fval))]);
            disp(['    wavelet space absolute err : ', num2str(err_wavelet)]);
            disp(['    wavelet space relative err : ', num2str(err_wavelet/max(abs(fval_analytic(:)))*100), ' %']);
        end
        
        %%
        % Check the realspace solution
        
        if num_dimensions <= 3
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
    
    %%
    % Plot results
    
    if mod(L,opts.plot_freq)==0 && ~opts.quiet
        
        figure(1000)
        
        if num_dimensions <= 3
            
            f_realspace_nD = singleD_to_multiD(num_dimensions,fval_realspace,nodes);
            if strcmp(opts.output_grid,'fixed')
                f_realspace_nD = ...
                    remove_duplicates(num_dimensions,f_realspace_nD,nodes_nodups,nodes_count);
            end

            f_realspace_analytic_nD = get_analytic_realspace_solution_D(pde,opts,coord_nodups,t+dt);
            
            element_coordinates = [];
            if opts.use_oldhash
            else
                element_coordinates = get_sparse_grid_coordinates(pde,opts,hash_table);
            end
            plot_fval(pde,nodes_nodups,f_realspace_nD,f_realspace_analytic_nD,element_coordinates);
           
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
    
    count=count+1;
    t1 = toc;
    if ~opts.quiet; disp(['Took ' num2str(t1) ' [s]']); end
    
    %%
    % Save output
    
    if opts.save_output && (mod(L,opts.save_freq)==0 || L==num_steps)
        [status, msg, msgID] = mkdir([root_directory,'/output']);
        if isempty(opts.output_filename_id)
            adapt_str = 'n';
            if opts.adapt; adapt_str=num2str(opts.adapt_threshold,'%1.1e'); end
            filename_str = ['-l',replace(num2str(opts.lev),' ','') ...
                '-d',num2str(opts.deg),'-',opts.grid_type,'-dt',num2str(dt),...
                '-adapt-',adapt_str];
        else           
            filename_str = opts.output_filename_id;           
        end
        fName = append(root_directory,"/output/asgard-out",filename_str,".mat");
 
        f_realspace_nD = singleD_to_multiD(num_dimensions,fval_realspace,nodes); 
        if strcmp(opts.output_grid,'fixed')
            f_realspace_nD = ...
                remove_duplicates(num_dimensions,f_realspace_nD,nodes_nodups,nodes_count);
        end
        
        if num_dimensions >= 1
            f_realspace_nD_t{L+1} = f_realspace_nD;
        else
            error('Save output for num_dimensions >3 not yet implemented');
        end
        fval_t{L+1} = fval;
        time_array(L+1) = t+dt;
        wall_clock_time(L+1) = toc;
       
        save(fName,'pde','opts','dt','f_realspace_nD_t','fval_t','nodes','time_array','hash_table','wall_clock_time');

    end
    
    t = t + dt;
    
end

end

