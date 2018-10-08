%function [f2d] = convertFK6DtoRealSpace()


isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
idebug = 1;

% Read in output from FK6D and convert it to real space

% Make Octave Compatible
addpath(genpath(pwd));


filename = 'fval-MWT.h5';

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
    Lev = double(5);
    dt = 7.716*10^-11;
    Lmin = 0.0;
    Lmax = 1.0;
    Vmin = -4500000.0;
    Vmax = 4500000.0;
    x1 = linspace(0,1,32);
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

pde.params.Deg = Deg;

[nW,nT] = size(fval_t);
t = (0:nT-1)*dt;

LevX = Lev;
LevV = Lev;

Dim = 2;

[HASH,HASHInv] = HashTable(Lev,Dim);
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

if (idebug >= 1),
    %steps = [0:1:1000]+1;
    steps = [0:1:(nT-1)/1]+1;
else
    steps = [0:1:(nT-1)/1]+1;
end;
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

[fx_loc] = EvalWavPoint4(Lmin,Lmax,Lev,Deg,xx);
[fv_loc] = EvalWavPoint4(Vmin,Vmax,Lev,Deg,vv);


doTransform = 1;
startFromLast = 0;

cnt = 1;

if doTransform
    more off;
    
    time_main = tic();
    
    f2d_t = zeros(numel(xx), numel(vv), numel(steps));
    
    nHash = numel(HASHInv),
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
        
        K1 = calc_overlap_index( i1, n1, LevV );
        K2 = calc_overlap_index( i2, n2, LevX );

        if (idebug >= 1),
          K1_find = find( sum(fv_loc(index_I1,:) ~= 0, 1 ));
          K2_find = find( sum(fx_loc(index_I2,:) ~= 0, 1 ));
          isok_K1 = (numel(K1) == numel(K1_find)) && all( K1_find == K1 );
          isok_K2 = (numel(K2) == numel(K2_find)) && all( K2_find == K2 );
          if (~isok_K1),
            disp(sprintf('n1=%d, i1=%d', n1,i1));
            for i=1:numel(K1),
              disp(sprintf('K1(%d) = %g', i,K1(i)));
            end;
            for i=1:numel(K1_find),
              disp(sprintf('K1_find(%d) = %g ',i,K1_find(i)));
            end;
          end;
          if (~isok_K2),
            disp(sprintf('n2=%d,i2=%d', n2,i2));
            for i=1:numel(K2),
              disp(sprintf('K2(%d) = %g', i,K2(i)));
            end;
            for i=1:numel(K2_find),
              disp(sprintf('K2_find(%d) = %g ',i,K2_find(i)));
            end;
          end;
         end;


        
        Amat = transpose(full(fv_loc(index_I1,K1)));
        Bmat = transpose(full(fx_loc(index_I2,K2)));
        Xmat = fval_t(  Index, 1:nsteps);
        
        % -----------------------------------------------
        % Ymat = Bmat * Xmat * transpose(Amat)
        % Ymat(K2,K1) = transpose(fx_loc(index_I2, K2)) *
        %                  Xmat(Index_I2,Index_I1) *
        %                    ( fv_loc(index_I1,K1) )
        % -----------------------------------------------
        
        Ymat = kron_mult2(Amat,Bmat, Xmat );
        
        f2d_t(K2, K1, 1:nsteps) = f2d_t(K2, K1, 1:nsteps) + ...
            reshape( Ymat, [numel(K2), numel(K1), nsteps]);
        
        if (idebug >= 1),
            disp(sprintf('i=%d,n1=%d,n2=%d, time=%g', i,n1,n2,toc(time_main)));
        end;
    end;
    
    elapsed_time = toc( time_main );
    disp(sprintf('elapsed time for %d time steps is %g sec', ...
        nsteps,   elapsed_time ));
    
    if (idebug >= 1),
        % double check
        f2d_t_cal = f2d_t;
        clear f2d_t;
        load('f2d_t.mat', 'f2d_t');
        abserr = abs( f2d_t_cal(:,:,1:nsteps) -  f2d_t(:,:,1:nsteps) ) ;
        relerr = abserr ./ max( abs(f2d_t_cal(:,:,1:nsteps)), abs(f2d_t(:,:,1:nsteps) ));
        max_abserr = max( abserr(:) );
        max_relerr = max( relerr(:) );
        disp(sprintf('f2d_t: max_abserr =%g, max_relerr=%g', ...
            max_abserr,     max_relerr ));
    end;
end;







%% Begin plotting
for cnt=1:nsteps,
    n_t(:,cnt) = sum(f2d_t(:,:,cnt),2) * dx * dv;
end;

F(nT) = struct('cdata',[],'colormap',[]);

df2d_t = f2d_t-f2d_t(:,:,1);
dn_t = n_t-n_t(:,1);

figure(1)
if (isOctave),
    contourf(dn_t,30);
    colormap(hot);
else
    colormap(redblue);
end;

dRange = max(abs(dn_t(:)));
caxis([-dRange,+dRange]);

c = hot;
c = flipud(c);

figure(2)
if (isOctave),
    max_tt = 1;
else
    max_tt = nS;
end;

for tt = 1:max_tt,
    ax1 = subplot(2,2,1);
    mesh(xx,vv,df2d_t(:,:,tt)','FaceColor','interp','EdgeColor','none');
    axis([Lmin Lmax Vmin Vmax])
    if (isOctave),
        %colormap(ax1,c);
    else
        %colormap(ax1,redblue);
    end;
    dRange = max(abs(df2d_t(:)));
    caxis([-dRange,+dRange]);
    colorbar;
    
    ax2 = subplot(2,2,2);
    mesh(xx,vv,f2d_t(:,:,tt)','FaceColor','interp','EdgeColor','none');
    %contourf(xx,vv,f2d_t(:,:,tt)','LineStyle','none');
    
    axis([Lmin Lmax Vmin Vmax])
    %colormap(ax2,c);
    %colorbar
    
    subplot(2,2,3);
    
    plot(x,dn_t(:,tt));
    maxn = max(dn_t(:));
    minn = min(dn_t(:));
    %axis([x1(1) x1(end) minn maxn]);
    
    title('Density Perturbation');
    
    subplot(2,2,4);
    plot(x1,E_t(:,tt));
    maxE = max(abs(E_t(:)))+1e-8;
    axis([x1(1) x1(end) -maxE maxE]);
    
    if (isOctave),
        %  ---------------------------------
        %  getframe not available in octave
        %  ---------------------------------
    else
        F(tt) = getframe(gcf);
    end;
    title ('Electric Field');
    
    pause(0.01)
    
    %F(tt) = getframe(gcf);
    
end


