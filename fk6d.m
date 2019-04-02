%% MATLAB (reference) version of the ASGarD solver

function [err,fval,fval_realspace] = fk6d(pde,lev,deg,TEND,quiet,compression,implicit,gridType,useConnectivity)

format short e
folder = fileparts(which(mfilename));
addpath(genpath(folder));

%% Load PDE and runtime defaults
runtimeDefaults

%% Shortcuts (some of this will go away soon)
% Named domain ranges
if nDims==2
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
pde.CFL = 0.1;
dt = pde.set_dt(pde);
if ~quiet; disp(sprintf('dt = %g', dt )); end

%% Construct the 1D multi-wavelet transform for each dimension.
for d=1:nDims
    pde.dimensions{d}.FMWT = OperatorTwoScale(pde.dimensions{d}.deg,2^pde.dimensions{d}.lev);
end

%% Construct the Element (Hash) tables.
if ~quiet; disp('Constructing hash and inverse hash tables'); end
%[HASH,HASHInv,index1D] = HashTable2(Lev,Dim); % I dont recall what the index1D information was added here for?
[HASH,HASHInv] = HashTable(lev,nDims,gridType);
nHash = numel(HASHInv);

%% Construct the connectivity.
if runTimeOpts.useConnectivity
    if ~quiet; disp('Constructing connectivity table'); end
    connectivity = ConnectnD(nDims,HASH,HASHInv,lev,lev);
else
    connectivity = [];
end

%% Generate initial conditions (both 1D and multi-D).
if ~quiet; disp('Calculate 2D initial condition on the sparse-grid'); end
fval = initial_condition_vector(HASHInv,pde,0);

%% Construct the time-independent coefficient matrices
% The original way
if ~quiet; disp('Calculate time independent matrix coefficients'); end
if nDims==2
    [vMassV,GradV,GradX,DeltaX,FluxX,FluxV] = matrix_coeff_TI(LevX,LevV,deg,Lmin,Lmax,Vmin,Vmax,...
        pde.dimensions{1}.FMWT,pde.dimensions{2}.FMWT);
end
%%
% The generalized PDE spec way
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
    A_data = GlobalMatrixSG_SlowVersion(pde,runTimeOpts,HASHInv,connectivity,deg);
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
for d=1:nDims
    [Meval{d},nodes{d}] = matrix_plot_D(pde.dimensions{d});
end
if nDims==2
    if ~quiet; disp('Plotting intial condition'); end
    [Meval_v,v_node,Meval_x,x_node]=matrix_plot(LevX,LevV,deg,Lmin,Lmax,Vmin,Vmax,...
        pde.dimensions{1}.FMWT,pde.dimensions{2}.FMWT);
    tol = 1e-15;
    assert(norm(full(Meval{1}-Meval_v))<tol);
    assert(norm(nodes{1}-v_node)<tol);
    assert(norm(full(Meval{2}-Meval_x))<tol);
    assert(norm(nodes{2}-x_node)<tol);
end

%%
% Construct a n-D coordinate array
% TODO : generalize to dimension better.

if nDims <= 3
    
    if nDims ==1
        [xx1] = ndgrid(nodes{1});
        coord = {xx1};
    end
    if nDims==2
        [xx1,xx2] = ndgrid(nodes{2},nodes{1});
        coord = {xx2,xx1};
    end
    if nDims==3
        [xx1,xx2,xx3] = ndgrid(nodes{3},nodes{2},nodes{1});
        coord = {xx3,xx2,xx1};
    end
    if nDims==6
        [xx1,xx2,xx3,xx4,xx5,xx6] = ndgrid(nodes{6},nodes{5},nodes{4},nodes{3},nodes{2},nodes{1});
        coord = {xx6,xx5,xx4,xx3,xx2,xx1};
    end
    
end

% %%
% % Try transforming a known 3D function to wavelet space and then back again.
%
% fa = getAnalyticSolution_D(coord,5*dt,pde);
% fa_wSpace = exact_solution_vector(HASHInv,pde,5*dt);
% fa_rSpace = reshape(Multi_2D_D(Meval,fa_wSpace,HASHInv,pde),size(fa));
%
% norm(fa(:)-fa_rSpace(:))/norm(fa(:))*100
%
% sy=5;sz=10;
% figure
% hold on
% plot(fa(:,sy,sz));
% plot(fa_rSpace(:,sy,sz));
%
% figure
% subplot(2,2,1)
% contour(fa_rSpace(:,:,sz));
% hold on
% subplot(2,2,2)
% contour(fa(:,:,sz));
% subplot(2,2,3)
% fa2=permute(fa,[3,1,2]);
% contour(fa2(:,:,sz));
% subplot(2,2,4)
% fa2=permute(fa,[3,2,1]);
% contour(fa2(:,:,sz));

%% Plot initial condition
if nDims <=3
    
    %%
    % Get the real space solution
    fval_realspace = Multi_2D_D(Meval,fval,HASHInv,pde);
    fval_realspace_analytic = getAnalyticSolution_D(coord,0,pde);

    if norm(fval_realspace) > 0
        plot_fval(pde,nodes,fval_realspace,fval_realspace_analytic);
    end
    
end

if nDims==2
    if ~quiet
        %%
        % Transform from wavelet space to real space
        tmp = Multi_2D(Meval_v,Meval_x,fval,HASHInv,lev,deg);
        figure(1000)
        
        f2d0 = reshape(tmp,deg*2^LevX,deg*2^LevV)';
        
        [xx,vv]=meshgrid(x_node,v_node);
        
        ax1 = subplot(1,2,1);
        mesh(xx,vv,f2d0,'FaceColor','interp','EdgeColor','none');
        axis([Lmin Lmax Vmin Vmax])
        %caxis([-range1 +range1]);
        title('df');
        
        ax2 = subplot(1,2,2);
        mesh(xx,vv,f2d0,'FaceColor','interp','EdgeColor','none');
        axis([Lmin Lmax Vmin Vmax])
        %caxis([range2n +range2]);
        title('f');
    end
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
    timeStr = sprintf('Step %i of %i',L,nsteps);
    
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
    
    
    %     ax3 = subplot(2,2,3);
    %     plot(x_node,Meval_x*u,'r-o')
    %     title(['time = ',num2str(dt*L)])
    %     ax4 = subplot(2,2,4);
    %     plot(x_node,Meval_x*E,'r-o')
    
    
    if ~quiet; disp('    Calculate time dependent matrix coeffs'); end
    if nDims==2
        if (pde.applySpecifiedE | pde.solvePoisson)
            
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
    
    %%% Update A_encode for time-dependent coefficient matricies.
    if ~quiet; disp('    Generate A_encode for time-dependent coeffs'); end
    if runTimeOpts.compression == 3
        % the new matrix construction is as _newCon, only works for
        % compression= 3
        %         B_encode = GlobalMatrixSG(GradV,EMassX,HASHInv,Con2D,Deg);
        B_encode=GlobalMatrixSG_newCon(GradV,EMassX,HASH,lev,deg);
        C_encode=[A_encode B_encode];
    else
        
    end
    
    %%% Advance Vlasov in time with RK3 time stepping method.
    if ~quiet; disp('    RK3 time step'); end
    if runTimeOpts.compression == 3
        fval = TimeAdvance(pde,runTimeOpts,C_encode,fval,time(count),dt,deg,HASHInv);
    else
        
        if nDims==2
            A_data.GradX     = pde.terms{1}{2}.coeff_mat;
            A_data.vMassV    = pde.terms{1}{1}.coeff_mat;
            A_data.EMassX    = pde.terms{2}{2}.coeff_mat;
            A_data.GradV     = pde.terms{2}{1}.coeff_mat;
            
            A_data.FluxX = FluxX;
            A_data.FluxV = FluxV;
        end
        
        % Write the A_data structure components for use in HPC version.
        write_A_data = 0;
        if write_A_data && L==1; write_A_data_to_file(A_data,lev,deg); end
        
        if nDims~=2
            Vmax = 0;
            Emax = 0; % These are only used in the global LF flux
        end
        fval = TimeAdvance(pde,runTimeOpts,A_data,fval,time(count),dt,deg,HASHInv,Vmax,Emax);
        
    end
    
    %%% Write the present fval to file.
    if write_fval; write_fval_to_file(fval,lev,deg,L); end
    
    %%% Write data for FK6D test
    
    %     fname = ['tests/vlasov4_time_5_3/fval_',num2str(L,'%3.3i'),'.dat'];
    %     fd = fopen(fname,'w'); % where file.dat is the name you want to save to
    %     fwrite(fd,full(fval),'double'); % where U is the vector/matrix you want to store, double is the typename
    %     fclose(fd);
    
    if nDims <=3
        
        %%
        % Get the real space solution
        fval_realspace = Multi_2D_D(Meval,fval,HASHInv,pde);
        
        %%
        % Try with function convertToRealSpace
        
        tryConvertToRealSpace = 1;
        if tryConvertToRealSpace
            LminB = zeros(1,nDims);
            LmaxB = zeros(1,nDims);
            for d=1:nDims
                LminB(d) = pde.dimensions{d}.domainMin;
                LmaxB(d) = pde.dimensions{d}.domainMax;
            end
            fval_realspaceB = converttoRealSpace(nDims,lev,deg,gridType,LminB,LmaxB,fval,lev);
        end
        
    end
    
    %%
    % Check against known solution
    if pde.checkAnalytic
        
        %%
        % Check the wavelet space solution with the analytic solution
        
        fval_analytic = exact_solution_vector(HASHInv,pde,L*dt);
        err_wavelet = sqrt(mean((fval(:) - fval_analytic(:)).^2));
        disp(['    wavelet space absolute err : ', num2str(err_wavelet)]);
        disp(['    wavelet space relative err : ', num2str(err_wavelet/max(abs(fval_analytic(:)))*100), ' %']);
        
        if nDims <= 3
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
        
        plot_fval(pde,nodes,fval_realspace,fval_realspace_analytic);
        
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
    
end

end

