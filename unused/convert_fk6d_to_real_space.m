% Run with no arguments to test.

function [f2d_t] = convertFK6DtoRealSpace(filename, gridType)
if ~exist('gridType','var') || isempty(gridType)
    gridType = 'SG'%'FG';
else
    if strcmp(gridType,'SG') || strcmp(gridType,'FG')
    else
        error("gridType must be set to 'SG' or 'FG'"); 
    end
end
addpath(genpath('./'));

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
idebug = 1;
use_find = 0;

if ~exist('filename','var') || isempty(filename)
    filename = 'tests/convertFK6DtoRealSpace/convert-test-fval-MWT.h5';
    filename2 = 'tests/convertFK6DtoRealSpace/convert-test-f2d_t.mat';
end

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

if (isOctave),
    % -----------------------------
    % perform octave test framework
    % -----------------------------
    data = load(filename, '-hdf5');
    fval_t = data.fval;
    E_t = data.E;
    % ------------------------------------
    % run Python script to grab attributes
    % ------------------------------------
    Deg = double(4);
    Lev = double(6);
    %dt = 7.716*10^-11;
    dt = 0.002797;
    Lmin = 0.0;
    Lmax = 20.944;
    Vmin = -13.0;
    Vmax = 13.0;
    x1 = linspace(0,1,2^Lev);
else
    timeStride = 1;
    info = h5info(filename);
    f_size = info.Datasets(2).Dataspace.Size;
    fval_t = h5read(filename, '/fval',[1 1],[f_size(1) f_size(2)/timeStride],[1 timeStride]);
    E_t = h5read(filename, '/E');
    
    Deg = double(h5readatt(filename, '/fval/', 'deg'));
    Lev = double(h5readatt(filename, '/fval/', 'lev'));
    dt = h5readatt(filename, '/fval/', 'dt');
    Lmin = h5readatt(filename, '/fval/', 'Lmin');
    Lmax = h5readatt(filename, '/fval/', 'Lmax');
    Vmin = h5readatt(filename, '/fval/', 'Vmin');
    Vmax = h5readatt(filename, '/fval/', 'Vmax');
    x1 = h5readatt(filename, '/fval/', 'x1');
end;

Lev_solution = Lev;
Lev_viz = Lev + 2;
Lev = Lev_viz;

if (idebug >= 1),
    disp(sprintf('Lev_solution=%d, Lev_viz=%d', ...
		  Lev_solution, Lev_viz));
end;

pde.params.Deg = Deg;

[nW,nT] = size(fval_t);
t = (0:nT-1)*dt;

LevX = Lev;
LevV = Lev;

Dim = 2;

[HASH,HASHInv] = HashTable(Lev_solution,Dim,gridType);
if (idebug >= 1),
    disp(sprintf('numel(HASH)=%d, numel(HASHInv)=%d', ...
        numel(HASH),    numel(HASHInv) ));
end;

n = 2^Lev*(Deg-1);

nX = n;
nY = n;

x = linspace(Lmin,Lmax,nX);
v = linspace(Vmin,Vmax,nY);

dx = x(2)-x(1);
dv = v(2)-v(1);


steps = [0:1:(nT-1)/1]+1;


nS = numel(steps);
nsteps = numel(steps);

n_t = zeros(nX,nS);

%% 2D Reverse MWT

% Evaluation starts from here
% Method :: evaluate every basis at the given points (xx) and (vv)
% fx_loc = phi_x (xx) with size DoF1DxNum_xx
% fv_loc = phi_y (yy) with size DoF1DxNum_vv
% Loop over every point with (xx_i,yy_j) and combine the value
%   for fx_loc(:,i) and fv_loc(:,j)

% if we have computed fval, we can just start with following
% Input:: Lstart, Lend, Lev, Deg, xx, yy

xx = linspace(Lmin,Lmax,nX)';
vv = linspace(Vmin,Vmax,nY)';

% ----------------------------------------------
% Note that fx_loc and fv_loc are small matrices
% may consider fully dense storage
% sparsity ratio about 18.75 percent
% ----------------------------------------------
[fx_loc] = EvalWavPoint4(Lmin,Lmax,Lev,Deg,xx);
[fv_loc] = EvalWavPoint4(Vmin,Vmax,Lev,Deg,vv);

doTransform = 1;
startFromLast = 0;

cnt = 1;

if doTransform
    more off;
    
    time_main = tic();
    
    f2d_t = zeros(numel(xx), numel(vv), numel(steps));
    
    nHash = numel(HASHInv);
    for i=1:nHash,
        ll = HASHInv{i};
        n1 = ll(1);
        n2 = ll(2);
        i1 = ll(3);
        i2 = ll(4);
        I1 = ll(5);
        I2 = ll(6);
        index_1 = I1;
        index_2 = I2;
        index_I1 = (I1-1)*Deg + (1:Deg);
        index_I2 = (I2-1)*Deg + (1:Deg);
        Index = (i-1)*(Deg*Deg) + (1:(Deg*Deg));
        
        isok = (I1 == LevCell2index(n1,i1)) & ...
            (I2 == LevCell2index(n2,i2));
        if (~isok),
            error(sprintf('i=%d,I1=%d,I2=%d, n1=%d,i1=%d, n2=%d,i2=%d', ...
                i,   I1,   I2,    n1,   i1,    n2,   i2));
        end;
        
        % --------------------------------------------------
        % Note K1 and K2 should have contiguous index values
        % --------------------------------------------------
        if (use_find),
            K1 = find( sum(fv_loc(index_I1,:) ~= 0, 1 ));
            K2 = find( sum(fx_loc(index_I2,:) ~= 0, 1 ));
        else
            [K1_start,K1_end] = calc_overlap_index(n1,i1, Lev, Deg);
            [K2_start,K2_end] = calc_overlap_index(n2,i2, Lev, Deg);
            K1 = K1_start:K1_end;
            K2 = K2_start:K2_end;
        end;
        
        % --------------------------------------------
        % note convert to full storage to use efficient
        % DGEMM operations
        % --------------------------------------------
        Amat = transpose(full(fv_loc(index_I1,K1)));
        Bmat = transpose(full(fx_loc(index_I2,K2)));
        Xmat = fval_t(  Index, 1:nsteps);
        
        % -----------------------------------------------
        % Ymat = Bmat * Xmat * transpose(Amat)
        % Ymat(K2,K1) = transpose(fx_loc(index_I2, K2)) *
        %                  Xmat(Index_I2,Index_I1) *
        %                    ( fv_loc(index_I1,K1) )
        % -----------------------------------------------
        
        Ymat = kronmult2(Amat,Bmat, Xmat );
        
        f2d_t(K2, K1, 1:nsteps) = f2d_t(K2, K1, 1:nsteps) + ...
            reshape( Ymat, [numel(K2), numel(K1), nsteps]);
        
        if (idebug >= 2),
            disp(sprintf('i=%d,n1=%d,n2=%d, size(K1)=%g, size(K2)=%g, time=%g', ...
                i,n1,n2,  numel(K1), numel(K2), toc(time_main)));
            
            % --------------------------------------------------
            % Note K1 and K2 should have contiguous index values
            % --------------------------------------------------
            max_K1 = max(K1(:)); min_K1 = min(K1(:));
            is_contiguous_K1 = numel(K1) == (max_K1-min_K1+1);
            max_K2 = max(K2(:)); min_K2 = min(K2(:));
            is_contiguous_K2 = numel(K2) == (max_K2-min_K2+1);
            
            [K1_start,K1_end] = calc_overlap_index(n1,i1, Lev, Deg);
            [K2_start,K2_end] = calc_overlap_index(n2,i2, Lev, Deg);
            
            isok = (K1_start == min_K1) && (K1_end == max_K1) && ...
                (K2_start == min_K2) && (K2_end == max_K2);
            
            
            
            if (~isok),
                disp(sprintf('** error ** n1=%d,i1=%d,K1=%d:%d,  K1_start:K1_end=%d:%d ', ...
                    n1,   i1,   min_K1,max_K1, K1_start,K1_end));
                disp(sprintf('** error ** n2=%d,i2=%d,K2=%d:%d,  K2_start:K2_end=%d:%d ', ...
                    n2,   i2,   min_K2,max_K2, K2_start,K2_end));
            end;
            
            
            if (~is_contiguous_K1),
                disp(sprintf('numel(K1)=%g  ** K1 is not contiguous ** ',...
                    numel(K1) ));
            end;
            
            if (~is_contiguous_K2),
                disp(sprintf('numel(K2)=%g  ** K2 is not contiguous ** ',...
                    numel(K2) ));
            end;
        end;
    end;
    
    elapsed_time = toc( time_main );
    disp(sprintf('elapsed time for %d time steps is %g sec', ...
        nsteps,   elapsed_time ));
    
    if (idebug >= 1),
      if (Lev_viz == Lev_solution),
        % double check
        f2d_t_cal = f2d_t;
        clear f2d_t;
        load(filename2, 'f2d_t');
        abserr = zeros(size(f2d_t,1),size(f2d_t,2));
        relerr = zeros(size(f2d_t,1),size(f2d_t,2));
        max_abserr = 0;
        max_relerr = 0;
        for istep=1:nsteps,
            abserr = abs( f2d_t_cal(:,:,istep) -  f2d_t(:,:,istep) ) ;
            dnorm = max( abs(f2d_t_cal(:,:,istep)), abs(f2d_t(:,:,istep)) );
            dnorm = max( 1e-6*ones(size(dnorm)), dnorm );
            relerr = abserr ./ dnorm;
            
            max_abserr = max( max_abserr, max( abserr(:)) );
            max_relerr = max( max_relerr, max( relerr(:)) );
        end;
        disp(sprintf('f2d_t: max_abserr =%g, max_relerr=%g', ...
            max_abserr,     max_relerr ));
     end;
    end;
end;


%% Begin plotting


df2d_t = f2d_t-f2d_t(:,:,1);

figure(2)
if (isOctave),
    max_tt = 1;
else
    max_tt = nS;
end;

range1 = max(abs(df2d_t(:)));
range2 = max(abs(f2d_t(:)));

range2n = -range2;
if abs(min(f2d_t(:)))<0.5*range2
    range2n = 0;
end

rangeFac1 = 0.3;
rangeFac2 = 0.5;

writeVid = 0;
if writeVid
    vid = VideoWriter('fk6d','MPEG-4');
    open(vid);
end

for tt = 1:max_tt,
    
    ax1 = subplot(1,2,1);
    mesh(xx,vv,df2d_t(:,:,tt)','FaceColor','interp','EdgeColor','none');
    axis([Lmin Lmax Vmin Vmax])
    caxis([-range1 +range1]*rangeFac1);
    title('df');
    
    ax2 = subplot(1,2,2);
    mesh(xx,vv,f2d_t(:,:,tt)','FaceColor','interp','EdgeColor','none');
    axis([Lmin Lmax Vmin Vmax])
    caxis([range2n +range2]*rangeFac2);
    title('f');
    
    if writeVid
        frame = getframe(gcf);
        writeVideo(vid,frame);
    else
        pause(0.5);
    end
    
end

if writeVid
    close(vid);
end

end

