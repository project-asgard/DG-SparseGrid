%% MATLAB (reference) version of the ASGarD solver

function [err,fval,fval_realspace] = asgard (pde, varargin)

format short e
folder = fileparts(which(mfilename));
addpath(genpath(folder));

%% Load PDE and runtime defaults
runtime_defaults

%% Check PDE
pde = check_pde(pde);

%% Set time step.
dt = pde.set_dt(pde,opts.CFL);
if ~opts.quiet; disp(sprintf('dt = %g', dt )); end

%% Construct the Element (Hash) table.
if ~opts.quiet; disp('Constructing hash and inverse hash tables'); end

if opts.use_oldhash
    [HASH,hash_table] = hash_table_nD(pde.lev_vec, opts.grid_type);
else
    [elements, elements_idx]    = element_table (pde.lev_vec, pde.max_lev, opts.grid_type);
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
fval = initial_condition_vector(pde, opts, hash_table, 0);

%% Construct the time-independent coefficient matrices
if ~opts.quiet; disp('Calculate time independent matrix coefficients'); end
t = 0;
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

%% Construct a n-D coordinate array
coord = get_realspace_coords(pde,nodes);

%% Plot initial condition
if num_dimensions <=3
    
    %%
    % Get the real space solution
    fval_realspace = wavelet_to_realspace(pde,opts,Meval,fval,hash_table);

    fval_realspace_analytic = get_analytic_realspace_solution_D(coord,0,pde);
    
    if norm(fval_realspace) > 0 && ~opts.quiet
        plot_fval(pde,nodes,fval_realspace,fval_realspace_analytic);
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

%% Time Loop
count=1;
plotFreq = 1;
err = 0;
if ~opts.quiet; disp('Advancing time ...'); end
for L = 1:num_steps
    
    tic;
    time(count) = (L-1)*dt;
    timeStr = sprintf('Step %i of %i at %f seconds',L,num_steps,time(count));
    
    if ~opts.quiet; disp(timeStr); end
    Emax = 0;
    
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
        E = E * pde.Et(time(count),params);
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
    
    t = time(count);
    TD = 1;
    pde = get_coeff_mats(pde,t,TD);
    
    %%
    % Advance in time
    
    if ~opts.quiet; disp('    RK3 time step'); end
    
    %%
    % Write the A_data structure components for use in HPC version.
    write_A_data = 0;
    if write_A_data && L==1; write_A_data_to_file(A_data,lev,deg); end
    
    fval = time_advance(pde,opts,A_data,fval,time(count),dt,pde.deg,hash_table,[],[]);
    
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
        
        fval_analytic = exact_solution_vector(pde,opts,hash_table,L*dt);
        err_wavelet = sqrt(mean((fval(:) - fval_analytic(:)).^2));
        if ~opts.quiet       
            disp(['    wavelet space absolute err : ', num2str(err_wavelet)]);
            disp(['    wavelet space relative err : ', num2str(err_wavelet/max(abs(fval_analytic(:)))*100), ' %']);
        end
        
        %%
        % Check the realspace solution
        
        if num_dimensions <= 3
            fval_realspace_analytic = get_analytic_realspace_solution_D(coord,L*dt,pde);
            err_real = sqrt(mean((fval_realspace(:) - fval_realspace_analytic(:)).^2));
            if ~opts.quiet         
                disp(['    real space absolute err : ', num2str(err_real)]);
                disp(['    real space relative err : ', num2str(err_real/max(abs(fval_realspace_analytic(:)))*100), ' %']);
            end
        end
        
        err = err_wavelet;
    end
    
    %%
    % Plot results
    
    if mod(L,plotFreq)==0 && ~opts.quiet
        
        figure(1000)
        
        if num_dimensions <= 3
            plot_fval(pde,nodes,fval_realspace,fval_realspace_analytic);
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
    
    %%
    % Apply adaptivity
    
    if opts.adapt
        if ~opts.quiet; disp('Adapt grid ...'); end
        [pde,fval,hash_table,A_data,Meval,nodes,coord] = adapt(pde,opts,fval,hash_table,nodes,fval_realspace);
    end
    
end

end

