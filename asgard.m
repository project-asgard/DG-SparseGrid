%% MATLAB (reference) version of the ASGarD solver

function [err,fval,fval_realspace] = asgard(pde,lev,deg,TEND,quiet,compression,implicit,gridType,useConnectivity,CFL)

format short e
folder = fileparts(which(mfilename));
addpath(genpath(folder));

%% Load PDE and runtime defaults
runtimeDefaults

%% Check PDE
pde = checkPDE(pde);

%% Check dimensions
pde = checkAllDimensions(pde);

%% Check terms
pde = checkTerms(pde);

%% Shortcuts (some of this will go away soon)
% Named domain ranges
if num_dimensions==2
    Lmin = pde.dimensions{2}.domainMin;
    Lmax = pde.dimensions{2}.domainMax;
    Vmin = pde.dimensions{1}.domainMin;
    Vmax = pde.dimensions{1}.domainMax;
    
    %%
    % Level information.
    LevX = pde.dimensions{2}.lev;
    LevV = pde.dimensions{1}.lev;
end

%%
% Things to be removed
deg = pde.dimensions{1}.deg;
lev = pde.dimensions{1}.lev;
DimX = 1;
params = pde.params;

%% Set time step.
dt = pde.set_dt(pde);
if ~quiet; disp(sprintf('dt = %g', dt )); end

%% Construct the Element (Hash) tables.
if ~quiet; disp('Constructing hash and inverse hash tables'); end

pde.useHash  = 1;
pde.do_adapt = 0;

[HASH,HASHInv] = HashTable(pde,lev,num_dimensions,gridType); % TODO : move this call inside the if below.

if pde.useHash
else
    [elements, elements_idx]    = element_table (pde,opts);
    pde.elements                = elements;
    pde.elementsIDX             = elements_idx; % only to get the same order as the hash table
end

%% Construct the 1D multi-wavelet transform for each dimension.
for d=1:num_dimensions
    pde.dimensions{d}.FMWT = OperatorTwoScale(pde,d,deg,pde.dimensions{d}.lev);
end

%% Construct the connectivity.
if opts.useConnectivity
    if ~quiet; disp('Constructing connectivity table'); end
    connectivity = ConnectnD(num_dimensions,HASH,HASHInv,lev,lev);
else
    connectivity = [];
end

%% Generate initial conditions (both 1D and multi-D).
if ~quiet; disp('Calculate 2D initial condition on the sparse-grid'); end
fval = initial_condition_vector(HASHInv,pde,0);

%% Construct the time-independent coefficient matrices
if ~quiet; disp('Calculate time independent matrix coefficients'); end
t = 0;
TD = 0;
pde = getCoeffMats(pde,t,TD);

%% Construct A_encode / A_data time independent data structures.
if ~quiet; disp('Generate A_encode data structure for time independent coefficients'); end
if compression == 3
    % the new matrix construction is as _newCon, only works for
    % compression= 3
    %     A_encode=GlobalMatrixSG(vMassV,GradX,HASHInv,Con2D,Deg);
    A_encode=GlobalMatrixSG_newCon(vMassV,GradX,HASH,lev,deg,gridType);
else
    % A_data is constructed only once per grid refinement, so can be done
    % on the host side.
    A_data = GlobalMatrixSG_SlowVersion(pde,opts,HASHInv,connectivity,deg);
end

%% Construct Poisson matrix
% Compute the global matrix for spatial variables "x" by
%
% Poisson Solver: A_Poisson (Hash, Dim_x,k,LevX,DeltaX)
%
% Input: Hash, Dim_x,k,LevX,DeltaX,or nu, eps, CurlCurlX Output: A_Poisson
% Another Idea is to solve Poisson Equation on the finest full grid

if ~quiet; disp('Construct matrix for Poisson solve'); end
if pde.solvePoisson
    if DimX>1
        % Construct DeltaX for DimX
    else
        A_Poisson = DeltaX;
    end
end

%% Construct RMWT (Reverse Multi Wavelet Transform) in 2D
% Get the wavelet -> realspace transform matrices and realspace node
% locations for each dimension.
for d=1:num_dimensions
    [Meval{d},nodes{d}] = matrix_plot_D(pde.dimensions{d});
end

%% Construct a n-D coordinate array
coord = get_realspace_coords(pde,nodes);

%% Plot initial condition
if num_dimensions <=3
    
    %%
    % Get the real space solution
    fval_realspace = Multi_2D_D(pde,Meval,fval,HASHInv);
    fval_realspace_analytic = getAnalyticSolution_D(coord,0,pde);

    if norm(fval_realspace) > 0
        plot_fval(pde,nodes,fval_realspace,fval_realspace_analytic);
    end
    
    if pde.useHash
    else
        coordinates = get_sparse_grid_coordinates(pde);
    end
%     fval_realspace_SG = real_space_solution_at_coordinates_irregular(pde,fval,coordinates);
    
end

%% Write the initial condition to file.
write_fval = 0;
if write_fval; write_fval_to_file(fval,lev,deg,0); end

count=1;
plotFreq = 1;
err = 0;

%% Time Loop

if ~quiet; disp('Advancing time ...'); end
nsteps = max(1,floor( TEND/dt));
for L = 1:nsteps,
    
    tic;
    time(count) = (L-1)*dt;
    timeStr = sprintf('Step %i of %i at %f seconds',L,nsteps,time(count));
    
    if ~quiet; disp(timeStr); end
    Emax = 0;
    
    if pde.solvePoisson
        %%% Solve Poisson to get E (from 1-rho=1-int f dv)
        if ~quiet; disp('    Solve poisson to get E'); end
        %[E,u] = PoissonSolve2(LevX,Deg,Lmax,fval,A_Poisson,FMWT_COMP_x,Vmax,index1D);
        [E,u] = PoissonSolve(LevX,deg,Lmax,fval,A_Poisson,FMWT_COMP_x,Vmax);
        Emax = max(abs(Meval{2}*E)); % TODO : this clearly is problem dependent
    end
    
    if pde.applySpecifiedE
        %%% Apply specified E
        if ~quiet; disp('    Apply specified E'); end
        E = forwardMWT(LevX,deg,Lmin,Lmax,pde.Ex,pde.params);
        E = E * pde.Et(time(count),params);
        Emax = max(abs(Meval{2}*E)); % TODO : this clearly is problem dependent
    end
        
    if ~quiet; disp('    Calculate time dependent matrix coeffs'); end
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
    pde = getCoeffMats(pde,t,TD);
    
    %%
    % Update A_encode for time-dependent coefficient matricies.
    if ~quiet; disp('    Generate A_encode for time-dependent coeffs'); end
    if opts.compression == 3
        % the new matrix construction is as _newCon, only works for
        % compression= 3
        %         B_encode = GlobalMatrixSG(GradV,EMassX,HASHInv,Con2D,Deg);
        B_encode=GlobalMatrixSG_newCon(GradV,EMassX,HASH,lev,deg);
        C_encode=[A_encode B_encode];
    else
        
    end
    
    %%
    % Advance in time
    
    if ~quiet; disp('    RK3 time step'); end
    if opts.compression == 3
        fval = TimeAdvance(pde,opts,C_encode,fval,time(count),dt,deg,HASHInv);
    else       
        % Write the A_data structure components for use in HPC version.
        write_A_data = 0;
        if write_A_data && L==1; write_A_data_to_file(A_data,lev,deg); end
        
        if num_dimensions~=2
            Vmax = 0;
            Emax = 0; % These are only used in the global LF flux
        end
        fval = TimeAdvance(pde,opts,A_data,fval,time(count),dt,deg,HASHInv,Vmax,Emax);
        
    end
    
    %%% Write the present fval to file.
    if write_fval; write_fval_to_file(fval,lev,deg,L); end
    
    %%% Write data for FK6D test
    
    %     fname = ['tests/vlasov4_time_5_3/fval_',num2str(L,'%3.3i'),'.dat'];
    %     fd = fopen(fname,'w'); % where file.dat is the name you want to save to
    %     fwrite(fd,full(fval),'double'); % where U is the vector/matrix you want to store, double is the typename
    %     fclose(fd);
    
    if num_dimensions <=3
        
        %%
        % Get the real space solution
        fval_realspace = Multi_2D_D(pde,Meval,fval,HASHInv);
        
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
            fval_realspaceB = converttoRealSpace(pde,num_dimensions,lev,deg,gridType,LminB,LmaxB,fval,lev);
%             fval_realspace = fval_realspaceB;
        end
        
        if pde.useHash
        else
            coordinates = get_sparse_grid_coordinates(pde);
        end
%         fval_realspace_SG = real_space_solution_at_coordinates(pde,fval,coordinates);
        
    end
    
    %%
    % Check against known solution
    if pde.checkAnalytic
        
        %%
        % Check the wavelet space solution with the analytic solution
        
        fval_analytic = exact_solution_vector(pde,HASHInv,L*dt);
        err_wavelet = sqrt(mean((fval(:) - fval_analytic(:)).^2));
        disp(['    wavelet space absolute err : ', num2str(err_wavelet)]);
        disp(['    wavelet space relative err : ', num2str(err_wavelet/max(abs(fval_analytic(:)))*100), ' %']);
        
        if num_dimensions <= 3
            %%
            % Check the realspace solution
            
            fval_realspace_analytic = getAnalyticSolution_D(coord,L*dt,pde);
            
            err_real = sqrt(mean((fval_realspace(:) - fval_realspace_analytic(:)).^2));
            disp(['    real space absolute err : ', num2str(err_real)]);
            disp(['    real space relative err : ', num2str(err_real/max(abs(fval_realspace_analytic(:)))*100), ' %']);
        end
        
        err = err_wavelet;
    end
    
    %%
    % Plot results
    
    if mod(L,plotFreq)==0 && ~quiet
        
        figure(1000)
        
        if num_dimensions <= 3
            plot_fval(pde,nodes,fval_realspace,fval_realspace_analytic);
        end
        
    end
    
    count=count+1;
    t1 = toc;
    disp(['Took ' num2str(t1) ' [s]']);
    
    %%% Save output
    saveOutput = 0;
    if saveOutput
        stat = mkdir('output');
        fName = ['output/f2d-' sprintf('%04.4d',L) '.mat'];
        f2d = reshape(fval_realspace,deg*2^LevX,deg*2^LevV)';
        save(fName,'f2d','fval');
    end
    
    %%
    % Apply adaptivity
    if ~quiet; disp('Adapt grid ...'); end
    if pde.do_adapt
        [pde,fval,A_data,Meval,nodes,coord] = adapt(pde,opts,fval,HASHInv,connectivity,nodes,fval_realspace);
    end
    
end

end

