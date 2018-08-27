function [f2d] = convertFK6DtoRealSpace()

% Read in output from FK6D and convert it to real space

filename = 'fval-MWT.h5';

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

pde = Vlasov7;

pde.params.Deg = Deg;

[nW,nT] = size(fval_t);
t = (0:nT-1)*dt;

LevX = Lev;
LevV = Lev;

Dim = 2;

[HASH,HASHInv] = HashTable(Lev,Dim);

n = 2^Lev;

nX = n;
nY = n;

x = linspace(Lmin,Lmax,nX);
n_t = zeros(nX,nT);

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

doTransform = 1;
if doTransform
    
    f2d_t = zeros(nX,nY,nT);
    
    for tt=1:nT
        
        fval = fval_t(:,tt);
        
        disp(tt);
        
        % loop over all the points on xx and yy
        for i = 1:length(xx)
            for j = 1:length(vv)
                
                fxList = {fx_loc(:,i)};
                fyList = {fy_loc(:,j)};
                ftList = {1};
                
                val = combine_dimensions_2(fxList,fyList,ftList,HASHInv,pde);
                
                f2d_t(i,j,tt) = val'*fval;
                
            end
        end
        
        % Density moment (hack)
        
        n_t(:,tt) = sum(f2d_t(:,:,tt),1);
        
    end
    
    save('f2d_t.mat','f2d_t', 'n_t', 'x', 'x1', 'E_t', 't');
     
else
    
    load('f2d_t.mat');
end


%% Begin plotting

F(nT) = struct('cdata',[],'colormap',[]);

df2d_t = f2d_t-f2d_t(:,:,1);
dn_t = n_t-n_t(:,1);

for tt=1:nT
    subplot(3,1,1);
    mesh(xx,vv,df2d_t(:,:,tt),'FaceColor','flat','EdgeColor','none');
    %mesh(xx,vv,f2d_t(:,:,tt),'FaceColor','interp','EdgeColor','none');
    
    axis([Lmin Lmax Vmin Vmax])
    
    subplot(3,1,2);
    
    plot(x,dn_t(:,tt));
    maxn = max(dn_t(:));
    minn = min(dn_t(:));
    axis([x1(1) x1(end) minn maxn]);
    
    subplot(3,1,3);
    plot(x1,E_t(:,tt));
    maxE = max(abs(E_t(:)));
    axis([x1(1) x1(end) -maxE maxE]);
    
    F(tt) = getframe(gcf);
    
end

view(0,90)
colorbar

end
