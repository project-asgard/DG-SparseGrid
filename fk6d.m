%% MATLAB (reference) version of the ASGarD solver

function [err,fval,fval_realspace] = fk6d(pde,Lev,Deg,TEND,quiet,compression,implicit,gridType)

format short e
folder = fileparts(which(mfilename));
addpath(genpath(folder));

%% Load PDE and runtime defaults
runtimeDefaults

%% Shortcuts (some of this will go away soon)
% Named domain ranges
Lmin = pde.dimensions{2}.domainMin;
Lmax = pde.dimensions{2}.domainMax;
Vmin = pde.dimensions{1}.domainMin;
Vmax = pde.dimensions{1}.domainMax;

%%
% Level information.
LevX = pde.dimensions{2}.lev;
LevV = pde.dimensions{1}.lev;

%%
% Things to be removed
Deg = pde.dimensions{1}.deg;
Lev = pde.dimensions{1}.lev;
DimX = 1;
params = pde.params;

%% Set time step.
CFL = 0.1;
dt = Lmax/2^LevX/Vmax/(2*Deg+1)*CFL;
if ~quiet; disp(sprintf('dt = %g', dt )); end

%% Construct the 1D multi-wavelet transform for each dimension.
for d=1:nDims
    pde.dimensions{d}.FMWT = OperatorTwoScale(pde.dimensions{d}.deg,2^pde.dimensions{d}.lev);
end

%% Construct the Element (Hash) tables.
if ~quiet; disp('Constructing hash and inverse hash tables'); end
%[HASH,HASHInv,index1D] = HashTable2(Lev,Dim); % I dont recall what the index1D information was added here for?
[HASH,HASHInv] = HashTable(Lev,nDims,gridType);
nHash = numel(HASHInv);

%% Construct the connectivity.
if ~quiet; disp('Constructing connectivity table'); end
Con2D = ConnectnD(nDims,HASH,HASHInv,Lev,Lev);

%% Generate initial conditions (both 1D and multi-D).
if ~quiet; disp('Calculate 2D initial condition on the sparse-grid'); end
% fval = initial_condition_vector(fx,fv,HASHInv,pde);
fval = initial_condition_vector(HASHInv,pde,0);

%% Construct the time-independent coefficient matrices
% The original way
if ~quiet; disp('Calculate time independent matrix coefficients'); end
[vMassV,GradV,GradX,DeltaX,FluxX,FluxV] = matrix_coeff_TI(LevX,LevV,Deg,Lmin,Lmax,Vmin,Vmax,...
     pde.dimensions{1}.FMWT,pde.dimensions{2}.FMWT);
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
    A_encode=GlobalMatrixSG_newCon(vMassV,GradX,HASH,Lev,Deg);
else
    % A_data is constructed only once per grid refinement, so can be done
    % on the host side.
    A_data = GlobalMatrixSG_SlowVersion(HASHInv,Con2D,Deg,compression);
end

%% Construct Poisson matrix
% Compute the global matrix for spatial variables "x" by
%
% Poisson Solver: A_Poisson (Hash, Dim_x,k,LevX,DeltaX)
%
% Input: Hash, Dim_x,k,LevX,DeltaX,or nu, eps, CurlCurlX Output: A_Poisson
% Another Idea is to solve Poisson Equation on the finest full grid

if ~quiet; disp('Construct matrix for Poisson solve'); end
if DimX>1
    % Construct DeltaX for DimX
else
    A_Poisson = DeltaX; 
end

%% Construct RMWT (Reverse Multi Wavelet Transform) in 2D
% Needs to be generalized to multi-D
if ~quiet; disp('Plotting intial condition'); end
[Meval_v,v_node,Meval_x,x_node]=matrix_plot(LevX,LevV,Deg,Lmin,Lmax,Vmin,Vmax,...
    pde.dimensions{1}.FMWT,pde.dimensions{2}.FMWT);
[xx,vv]=meshgrid(x_node,v_node);

%% Plot initial condition
if ~quiet
    %%
    % Transform from wavelet space to real space
    tmp = Multi_2D(Meval_v,Meval_x,fval,HASHInv,Lev,Deg);
    figure(1000)
    
    f2d0 = reshape(tmp,Deg*2^LevX,Deg*2^LevV)';
    
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

%% Write the initial condition to file.
write_fval = 0;
if write_fval; write_fval_to_file(fval,Lev,Deg,0); end

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
    
    if pde.solvePoisson
        %%% Solve Poisson to get E (from 1-rho=1-int f dv)
        if ~quiet; disp('    Solve poisson to get E'); end
        %[E,u] = PoissonSolve2(LevX,Deg,Lmax,fval,A_Poisson,FMWT_COMP_x,Vmax,index1D);
        [E,u] = PoissonSolve(LevX,Deg,Lmax,fval,A_Poisson,FMWT_COMP_x,Vmax);
    end
    
    if pde.applySpecifiedE
        %%% Apply specified E
        if ~quiet; disp('    Apply specified E'); end
        E = forwardMWT(LevX,Deg,Lmin,Lmax,pde.Ex,pde.params);
        E = E * pde.Et(time(count),params);
    end
    
    Emax = max(abs(Meval_x*E)); % max value on each point for E
    
    %     ax3 = subplot(2,2,3);
    %     plot(x_node,Meval_x*u,'r-o')
    %     title(['time = ',num2str(dt*L)])
    %     ax4 = subplot(2,2,4);
    %     plot(x_node,Meval_x*E,'r-o')
    
    %%% Generate EMassX time dependent coefficient matrix.
    if ~quiet; disp('    Calculate time dependent matrix coeffs'); end
    EMassX = matrix_coeff_TD(LevX,Deg,Lmin,Lmax,E,pde.dimensions{1}.FMWT);
    
    %%
    % Set the dat portion of the EMassX part of E.d_dv term.
    
    pde.terms{2}{2}.dat = E;
    
    %%
    % Now construct the TD coeff_mats.
    
    t = time(count);
    TD = 1;
    pde = getCoeffMats(pde,t,TD);
    
    %% Test new PDE spec based generation of the coeff_matrices
    
    disp( [ 'GradX error : '  num2str(norm(pde.terms{1}{2}.coeff_mat - GradX)/norm(GradX)) ]);
    disp( [ 'vMassV error : ' num2str(norm(pde.terms{1}{1}.coeff_mat - vMassV)/norm(vMassV)) ]);
    disp( [ 'EMassX error : ' num2str(norm(pde.terms{2}{2}.coeff_mat - EMassX)/norm(EMassX)) ]);
    disp( [ 'GradV error : '  num2str(norm(pde.terms{2}{1}.coeff_mat - GradV)/norm(GradV)) ]);
    
    %%% Update A_encode for time-dependent coefficient matricies.
    if ~quiet; disp('    Generate A_encode for time-dependent coeffs'); end
    if compression == 3
    % the new matrix construction is as _newCon, only works for 
    % compression= 3
%         B_encode = GlobalMatrixSG(GradV,EMassX,HASHInv,Con2D,Deg);
        B_encode=GlobalMatrixSG_newCon(GradV,EMassX,HASH,Lev,Deg);
        C_encode=[A_encode B_encode];
    else
        
    end
    
    %%% Advance Vlasov in time with RK3 time stepping method.
    if ~quiet; disp('    RK3 time step'); end
    if compression == 3
        fval = TimeAdvance(C_encode,fval,time(count),dt,compression,Deg,pde,HASHInv);
    else
        
        A_data.GradX     = pde.terms{1}{1}.coeff_mat;
        A_data.vMassV    = pde.terms{1}{2}.coeff_mat;
        A_data.EMassX    = pde.terms{2}{1}.coeff_mat;
        A_data.GradV     = pde.terms{2}{2}.coeff_mat;
        
        %         A_data.vMassV    = vMassV;
        %         A_data.GradX     = GradX;
        %         A_data.GradV     = GradV;
        %         A_data.EMassX    = EMassX;
        
        A_data.FluxX = FluxX;
        A_data.FluxV = FluxV;
        
        % Write the A_data structure components for use in HPC version.
        write_A_data = 0;
        if write_A_data && L==1; write_A_data_to_file(A_data,Lev,Deg); end
        
        fval = TimeAdvance(A_data,fval,time(count),dt,compression,Deg,pde,HASHInv,Vmax,Emax);
        
    end
    
    %%% Write the present fval to file.
    if write_fval; write_fval_to_file(fval,Lev,Deg,L); end
    
    %%% Write data for FK6D test
    
    %     fname = ['tests/vlasov4_time_5_3/fval_',num2str(L,'%3.3i'),'.dat'];
    %     fd = fopen(fname,'w'); % where file.dat is the name you want to save to
    %     fwrite(fd,full(fval),'double'); % where U is the vector/matrix you want to store, double is the typename
    %     fclose(fd);
    
    %%% Plot results
    if mod(L,plotFreq)==0 && ~quiet
        
        figure(1000)
        
        tmp=Multi_2D(Meval_v,Meval_x,fval,HASHInv,Lev,Deg);
        
        f2d = reshape(tmp,Deg*2^LevX,Deg*2^LevV)';
        
        %         ax1 = subplot(1,2,1);
        ax1 = subplot(2,2,1);
        mesh(xx,vv,f2d-f2d0,'FaceColor','interp','EdgeColor','none');
        axis([Lmin Lmax Vmin Vmax])
        %caxis([-range1 +range1]);
        title('df');
        
        %         ax2 = subplot(1,2,2);
        ax2 = subplot(2,2,2);
        mesh(xx,vv,f2d,'FaceColor','interp','EdgeColor','none');
        axis([Lmin Lmax Vmin Vmax])
        %caxis([range2n +range2]);
        title('f');
        
        title(['f @ ', timeStr])
        pause (0.01)
    end
    
    %%% Get the real space solution
    fval_realspace = Multi_2D(Meval_v,Meval_x,fval,HASHInv,Lev,Deg);
    
    %%% Check against known solution
    if pde.checkAnalytic
        
        % Check the wavelet space solution with the analytic solution
        fval_analytic = exact_solution_vector(HASHInv,pde,L*dt);
        err_wavelet = sqrt(mean((fval(:) - fval_analytic(:)).^2));
        disp(['    wavelet space absolute err : ', num2str(err_wavelet)]);
        disp(['    wavelet space relative err : ', num2str(err_wavelet/max(abs(fval_analytic(:)))*100), ' %']);
        
        % Check the real space solution with the analytic solution
        f2d = reshape(fval_realspace,Deg*2^LevX,Deg*2^LevV)';
        f2d_analytic = pde.analytic_solution(xx,vv,L*dt);
        err_real = sqrt(mean((f2d(:) - f2d_analytic(:)).^2));
        disp(['    real space absolute err : ', num2str(err_real)]);
        disp(['    real space relative err : ', num2str(err_real/max(abs(f2d_analytic(:)))*100), ' %']);
        
        err = err_wavelet;
    end
    
    count=count+1;
    t1 = toc;
    disp(['Took ' num2str(t1) ' [s]']);
    
    %%% Save output
    saveOutput = 0;
    if saveOutput
        stat = mkdir('output');
        fName = ['output/f2d-' sprintf('%04.4d',L) '.mat'];
        f2d = reshape(fval_realspace,Deg*2^LevX,Deg*2^LevV)';
        save(fName,'f2d','fval');
    end
    
end

end

