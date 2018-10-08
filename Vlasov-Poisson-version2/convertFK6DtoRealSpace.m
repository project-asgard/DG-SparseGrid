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
    dt = 7.716*10^2-11;
    Lmin = 0.0;
    Lmax = 1.0;
    Vmin = -4500000.0;
    Vmax = 4500000.0;
    x1 = linspace(0,1,32);
  else
    fval_t = h5read(filename, '/fval');
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
  steps = [0:1:5]+1;
else
  steps = [0:1:(nT-1)/1]+1;
end;
nS = numel(steps);

n_t = zeros(nX,nS);

%% 2D Reverse MWT

% Evaluation starts from here
% Method :: evaluate every basis at the given points (xx) and (vv)
% fx_loc = phi_x (xx) with size DoF1DxNum_xx
% fy_loc = phi_y (yy) with size DoF1DxNum_vv
% Loop over every point with (xx_i,yy_j) and combine the value
%   for fx_loc(:,i) and fy_loc(:,j)

% if we have computed fval, we can just start with following
% Input:: Lstart, Lend, Lev, Deg, xx, yy

xx = linspace(Lmin,Lmax,nX)';
vv = linspace(Vmin,Vmax,nY)';

[fx_loc] = EvalWavPoint4(Lmin,Lmax,Lev,Deg,xx);
[fy_loc] = EvalWavPoint4(Vmin,Vmax,Lev,Deg,vv);


preFileName = ['valPre-',num2str(Deg),'-',num2str(Lev),'.mat'];

if exist(preFileName) == 2
     
    load(preFileName);
    
else
    % loop over all the points on xx and yy
    disp('PreCalculate basis functions');
    bbb=1;
    for i = 1:length(xx)
        disp([num2str(i),' of ', num2str(length(xx))]);
        for j = 1:length(vv)
            
            fxList = {fx_loc(:,i)};
            fyList = {fy_loc(:,j)};
            ftList = {1};
            
            valPre(:,bbb) = combine_dimensions_2(fxList,fyList,ftList,HASHInv,pde);
            
            %f2d_t(i,j,cnt) = val'*fval;
            
            bbb = bbb + 1;
            
        end
    end
    
    save(preFileName,'valPre');
    
end


doTransform = 1;
startFromLast = 0;

cnt = 1;

if doTransform
    more off;
    
    if startFromLast
        load('f2d_t.mat');
        disp(['Restarting from last : ', num2str(cnt)]);
    else
        f2d_t = zeros(numel(xx),numel(vv), numel(steps));
    end
    
    for tt = steps
        
        fval = fval_t(:,tt);
        
        disp(cnt);
        
        % loop over all the points on xx and yy
        bbb = 1;
        for i = 1:length(xx)
            for j = 1:length(vv)
                
                %                 fxList = {fx_loc(:,i)};
                %                 fyList = {fy_loc(:,j)};
                %                 ftList = {1};
                %
                %                 val(:,bbb) = combine_dimensions_2(fxList,fyList,ftList,HASHInv,pde);
                
                val = valPre(:,bbb);
                
                f2d_t(i,j,cnt) = val'*fval;
                
                bbb = bbb +1;
                
            end
        end
        
        % Density moment (hack)
        
        n_t(:,cnt) = sum(f2d_t(:,:,cnt),2) * dx * dv;
        
        cnt = cnt + 1;
        
        %save('f2d_t.mat','f2d_t', 'n_t', 'x', 'x1', 'E_t', 't', 'cnt');
        
    end
    
    %% save('f2d_t.mat','f2d_t', 'n_t', 'x', 'x1', 'E_t', 't', 'cnt');
    
else
    
    load('f2d_t.mat');
end


%% Begin plotting

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
    mesh(xx,vv,df2d_t(:,:,tt)','FaceColor','flat','EdgeColor','none');
    axis([Lmin Lmax Vmin Vmax])
    if (isOctave),
      colormap(ax1,c);
    else
      colormap(ax1,redblue);
    end;
    dRange = max(abs(df2d_t(:)));
    caxis([-dRange,+dRange]);
    colorbar;
    
    ax2 = subplot(2,2,2);
    mesh(xx,vv,f2d_t(:,:,tt)','FaceColor','flat','EdgeColor','none');

    axis([Lmin Lmax Vmin Vmax])
    colormap(ax2,c);
    colorbar
    
    subplot(2,2,3);
    
    plot(x,dn_t(:,tt));
    maxn = max(dn_t(:));
    minn = min(dn_t(:));
    axis([x1(1) x1(end) minn maxn]);
    
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



%end
