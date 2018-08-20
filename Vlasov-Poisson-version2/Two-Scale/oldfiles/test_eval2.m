% Test
clear
% close all
% clc
format long e

Deg = 3;
Lstart = -1;
Lend = 2;
Lev = 2;
Lmax = Lend-Lstart;
hx = Lmax/2^Lev;

pde.params.Deg = Deg;
params.A = 0;


func_x = @(x,params)(sin(pi*x));
func_y = @(y,params)(sin(pi*y));



[f_coef_x] = forwardMWT(Lev,Deg,Lstart,Lend,func_x,params);
[f_coef_y] = forwardMWT(Lev,Deg,Lstart,Lend,func_y,params);


Dim = 2;
[HASH,HASHInv] = HashTable(Lev,Dim);
nHash = numel(HASHInv);

% base on the HashTable to generate 2D coefficient vector fval
fxList = {f_coef_x};
fyList = {f_coef_y};
ftList = {1};
[fval] = combine_dimensions_2(fxList,fyList,ftList,HASHInv,pde);

% Evaluation starts from here
% Method 1-- evaluate every basis at the given points (xx) and (vv)
% fx_loc = phi_x (xx) with size DoF1DxNum_xx
% fy_loc = phi_y (yy) with size DoF1DxNum_vv
% Loop over every point with (xx_i,yy_j) and combine the value
%   for fx_loc(:,i) and fy_loc(:,j)
% xx = [Lstart:hx:Lend]';
% yy = [Lstart:hx:Lend]';

xx = [Lstart:0.05:Lend]';
yy = [Lstart:0.05:Lend]';

[fx_loc,fx_val,ix,aijx] = EvalWavPoint3(Lstart,Lend,Lev,Deg,f_coef_x,xx);
[fy_loc,fy_val,iy,aijy] = EvalWavPoint3(Lstart,Lend,Lev,Deg,f_coef_y,yy);

% % % check with 1D evaluation
% % figure(3);
% % plot(xx,fx_val,'r-o',xx,func_x(xx),'b-<')
% % max(abs(fx_val-func_x(xx)))
% % % figure;
% % % plot(xx,fval-func(xx),'r-o')
% % % [norm(fval-func(xx)) max(abs(fval-func(xx)))]
% % 
% % return

% loop over all the points on xx and yy
for i = 1:length(xx)
    for j = 1:length(yy)
        
        fxList = {fx_loc(:,i)};
        fyList = {fy_loc(:,j)};
        ftList = {1};
        
        val = combine_dimensions_2(fxList,fyList,ftList,HASHInv,pde);
        
        f2D(i,j) = val'*fval;
        
    end
end
['Done with Method 1']

% Method 2
% Loop over each basis and evaluate the value at every point
n_points = length(xx)*length(yy);
npx = length(xx);
npy = length(yy);
f2D2 = sparse(Deg^2*nHash,n_points);

for i = 1 : nHash
    ll=HASHInv{i};
    
    % 1D indices for (Lev1,Cell1)-->Index1,(Lev2,Cell2)-->Index2
    I1=ll(5);
    I2=ll(6);
    
    Index1 = Deg*(I1-1)+[1:Deg];
    Index2 = Deg*(I2-1)+[1:Deg];
    
    tmp_ix = {ix{Index1}};
    tmp_aijx = {aijx{Index1}};
    
    tmp_iy = {iy{Index2}};
    tmp_aijy = {aijy{Index2}};
    
    for k1 = 1:Deg
        for k2 = 1:Deg
            tmp2D = tmp_aijx{k1}'*tmp_aijy{k2};
            tmp2D = tmp2D(:);
            
            A_loc = sparse(tmp_ix{k1},ones(length(tmp_ix{k1}),1),tmp_aijx{k1},npx,1);
            B_loc = sparse(tmp_iy{k2},ones(length(tmp_iy{k2}),1),tmp_aijy{k2},npy,1);
            
            tmp_loc = kron(A_loc,B_loc);
            
            Index0 = Deg^2*(i-1)+Deg*(k1-1)+k2;
            
            ii = Index0*ones(1,npx*npy);
            jj = 1:length(xx)*length(yy);
            
            f2D2 = f2D2 + sparse(ii,jj,tmp_loc,Deg^2*nHash,n_points);
            
        end
        
    end

end
['Done with Method 2']

% loop over Hash
[XX,YY] = meshgrid(xx,yy);
nx = length(xx);ny = length(yy);

figure(1);
subplot(1,3,1)
surf(XX,YY,f2D);%reshape(fval,nx,ny))
subplot(1,3,2)
surf(XX,YY,sin(pi*XX).*sin(pi*YY));

uu = sin(pi*XX).*sin(pi*YY);

[norm(f2D(:)-uu(:)) max(abs(f2D(:)-uu(:)))]

fnew = fval'*f2D2;
figure(1)
subplot(1,3,3)
surf(XX,YY,reshape(fnew,nx,ny))

figure(2)
plot(reshape(fnew,nx,ny)-f2D)

% another way to do inverse wavelet transform

return


% for 2D sparse grids solution

