%% MATLAB (reference) version of the ASGarD solver

function [err,fval,fval_realspace] = asgard (pde, varargin)

format short e
folder = fileparts(which(mfilename));
addpath(genpath(folder));

%% Load PDE and runtime defaults
runtime_defaults

%% Reset any persistent variables
if opts.time_independent_A
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
    [HASH,hash_table] = hash_table_nD(pde.lev_vec, opts.grid_type);
else
    [elements, elements_idx]    = hash_table_sparse_nD (pde.lev_vec, pde.max_lev, opts.grid_type);
    hash_table.elements         = elements;
    hash_table.elements_idx     = elements_idx; % only to get the same order as the old hash table
end

%% Construct the 1D multi-wavelet transform for each dimension.
for d=1:num_dimensions
    pde.dimensions{d}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{d}.lev);
end

%% (Do not) Construct the connectivity.
connectivity = [];

%% Generate initial conditions (both 1D and multi-D).
if ~opts.quiet; disp('Calculate 2D initial condition on the sparse-grid'); end
t = 0;
fval = initial_condition_vector(pde, opts, hash_table, t);

%% Construct the time-independent coefficient matrices
if ~opts.quiet; disp('Calculate time independent matrix coefficients'); end
TD = 0;
pde = get_coeff_mats(pde,t,TD);

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
    [Meval{d},nodes{d}] = matrix_plot_D(pde,pde.dimensions{d});
end

%%%%%%%%%
save('Meval','Meval')
save('nodes','nodes')
%%%%%%%%%


%% Construct a n-D coordinate array
coord = get_realspace_coords(pde,nodes);

%% Plot initial condition
if num_dimensions <=3
    
    %%
    % Get the real space solution
    fval_realspace = wavelet_to_realspace(pde,opts,Meval,fval,hash_table);
    fval_realspace_analytic = get_analytic_realspace_solution_D(pde,opts,coord,t);
    
    if norm(fval_realspace) > 0 && ~opts.quiet
        plot_fval(pde,nodes,fval_realspace,fval_realspace_analytic,Meval);
    end
    
    if opts.use_oldhash
    else
        coordinates = get_sparse_grid_coordinates(pde, opts, hash_table);
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
            [pde,fval_tmp,hash_table,A_data,Meval,nodes,coord] ...
                = adapt(pde,opts,fval,hash_table,nodes, ...
                fval_realspace,0,1);
            if num_pre_adapt == numel(fval_tmp)
                keep_adapting_initial_condition = false;
            end
            clear fval_tmp;
            % reproject initial condition onto refined basis
            fval = initial_condition_vector(pde, opts, hash_table, t);
        end
        
        % coarsen
        [pde,~,hash_table,A_data,Meval,nodes,coord] ...
            = adapt(pde,opts,fval,hash_table,nodes, ...
            fval_realspace,1,0);
        % reproject onto coarsend basis
        fval = initial_condition_vector(pde, opts, hash_table, t);

    else
        if ~opts.quiet; disp('Checking if initial grid is sufficient for requested accuracy ...'); end
        
        % Check to ensure refinement is not required to start
        pre_refinement_num_DOF = length(fval);
        [~,fval_check] ...
            = adapt(pde,opts,fval,hash_table,nodes,fval_realspace,0,1);
        if (length(fval_check)>pre_refinement_num_DOF)
%             error('Initial grid was insifficient for requested accuracy');
        end
        clear fval_check;
    end
else
   if ~opts.quiet; disp(['Number of DOF : ', num2str(numel(fval))]); end
end

fval_analytic = exact_solution_vector(pde,opts,hash_table,t);
fval_realspace = wavelet_to_realspace(pde,opts,Meval,fval,hash_table);
fval_realspace_analytic = get_analytic_realspace_solution_D(pde,opts,coord,t);

err_wavelet = sqrt(mean((fval(:) - fval_analytic(:)).^2));
err_realspace = sqrt(mean((fval_realspace(:) - fval_realspace_analytic(:)).^2));
if ~opts.quiet
    disp(['    num_dof : ', num2str(numel(fval))]);
    disp(['    wavelet space absolute err : ', num2str(err_wavelet)]);
    disp(['    wavelet space relative err : ', num2str(err_wavelet/max(abs(fval_analytic(:)))*100), ' %']);
    disp(['    wavelet space absolute err (2-norm) : ', num2str(norm(fval-fval_analytic))]);
    disp(['    real space absolute err : ', num2str(err_realspace)]);
    disp(['    real space relative err : ', num2str(err_realspace/max(abs(fval_realspace_analytic(:)))*100), ' %']);
    disp(['    real space absolute err (2-norm) : ', num2str(norm(fval_realspace(:)-fval_realspace_analytic(:)))]);
end

%% Time Loop
count=1;
plotFreq = 1;
err = 1e9;
if ~opts.quiet; disp('Advancing time ...'); end
for L = 1:num_steps
    
    tic;
    timeStr = sprintf('Step %i of %i at %f seconds',L,num_steps,t);
    
    if ~opts.quiet; disp(timeStr); end
    Emax = 0;
    
    % Coarsen Grid
    if opts.adapt
        [pde,fval,hash_table,A_data,Meval,nodes,coord] ...
            = adapt(pde,opts,fval,hash_table,nodes,fval_realspace,1,0);
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
                
                EMassX = matrix_coeff_TD(LevX,deg,Lmin,Lmax,E,pde.dimensions{1}.FMWT);
                
                %%
                % Set the dat portion of the EMassX part of E.d_dv term.
                
                pde.terms{2}{2}.dat = E;
                
            end
        end
        
        
        %%
        % Now construct the TD coeff_mats.
        
        TD = 1;
        pde = get_coeff_mats(pde,t,TD);
        
        %%
        % Advance in time
        
        if ~opts.quiet; disp(['    Time step (dt=',num2str(dt),')']); end
        
        %%
        % Write the A_data structure components for use in HPC version.
        write_A_data = 0;
        if write_A_data && L==1; write_A_data_to_file(A_data,lev,deg); end
        
        fval_unstepped = fval;
        fval = time_advance(pde,opts,A_data,fval,t,dt,pde.deg,hash_table,[],[]);
        fval_realspace = wavelet_to_realspace(pde,opts,Meval,fval,hash_table);
        
        %%
        % Refine Grid - determine which elements to add, but reset fval to
        % fval_previous with those new elements added and time_advance
        % again
        
        if opts.adapt
            if ~opts.quiet; disp('Adapt grid ...'); end
            
            num_elements_0 = numel(fval);
            
            [pde,~,hash_table,A_data,Meval,nodes,coord,fval_unstepped_adapted] ...
                = adapt(pde,opts,fval,hash_table,nodes, ...
                fval_realspace,0,1,fval_unstepped);
            
            num_elements_adapted = numel(fval_unstepped_adapted);
            
            if num_elements_0 == num_elements_adapted
                if ~opts.quiet
                    disp(['No more adaption needed - advancing time']);
                    disp(['    t = ', num2str(t)]);
                end
                needs_adapting = false;
                fval = time_advance(pde,opts,A_data,fval_unstepped_adapted,t,dt,pde.deg,hash_table,[],[]);
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
    
    %%% Write data for FK6D test
    
    %     fname = ['tests/vlasov4_time_5_3/fval_',num2str(L,'%3.3i'),'.dat'];
    %     fd = fopen(fname,'w'); % where file.dat is the name you want to save to
    %     fwrite(fd,full(fval),'double'); % where U is the vector/matrix you want to store, double is the typename
    %     fclose(fd);
    
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
            %             fval_realspace = fval_realspaceB;
        end
        
        if opts.use_oldhash
        else
            coordinates = get_sparse_grid_coordinates(pde,opts,hash_table);
        end
        %         fval_realspace_SG = real_space_solution_at_coordinates(pde,fval,coordinates);
        
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
            fval_realspace_analytic = get_analytic_realspace_solution_D(pde,opts,coord,t+dt);
            err_real = sqrt(mean((fval_realspace(:) - fval_realspace_analytic(:)).^2));
            if ~opts.quiet         
                disp(['    real space absolute err : ', num2str(err_real)]);
                disp(['    real space relative err : ', num2str(err_real/max(abs(fval_realspace_analytic(:)))*100), ' %']);
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
    
    if mod(L,plotFreq)==0 && ~opts.quiet
        
        figure(1000)
        
        if num_dimensions <= 3
            plot_fval(pde,nodes,fval_realspace,fval_realspace_analytic,Meval,coordinates);
        end
        
    end
    
    count=count+1;
    t1 = toc;
    if ~opts.quiet; disp(['Took ' num2str(t1) ' [s]']); end
    
    %%
    % Save output
    
    saveOutput = 0;
    if saveOutput
        stat = mkdir('output');
        fName = ['output/f2d-' sprintf('%04.4d',L) '.mat'];
        f2d = reshape(fval_realspace,deg*2^LevX,deg*2^LevV)';
        save(fName,'f2d','fval');
    end
    
    t = t + dt;
    
end

end

