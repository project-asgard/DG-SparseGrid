% Test
clear
close all
% clc
format long e

Deg = 2;
Lstart = -1;
Lend = 2;
Lev = 7;
Lmax = 1;


func_x = @(x,params)(sin(pi*x));
func_y = @(y,params)(sin(pi*y));

params.A = 0;
[f_coef_x] = forwardMWT(Lev,Deg,Lstart,Lend,func_x,params);
[f_coef_y] = forwardMWT(Lev,Deg,Lstart,Lend,func_y,params);

xx = [Lstart:0.1:Lend]';
yy = [Lstart:0.1:Lend]';
[fx_loc,fx_val] = EvalWavPoint3(Lstart,Lend,Lev,Deg,f_coef_x,xx);
[fy_loc,fy_val] = EvalWavPoint3(Lstart,Lend,Lev,Deg,f_coef_y,yy);

Dim = 2;
[HASH,HASHInv] = HashTable(Lev,Dim);
nHash = numel(HASHInv);


pde.params.Deg = Deg;

fxList = {f_coef_x};
fyList = {f_coef_y};
ftList = {1};
[fval] = combine_dimensions_2(fxList,fyList,ftList,HASHInv,pde);

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

% loop over Hash
[XX,YY] = meshgrid(xx,yy);
nx = length(xx);ny = length(yy);

figure;
subplot(1,2,1)
surf(XX,YY,f2D);%reshape(fval,nx,ny))
subplot(1,2,2)
surf(XX,YY,sin(pi*XX).*sin(pi*YY));

uu = sin(pi*XX).*sin(pi*YY);

[norm(f2D(:)-uu(:)) max(abs(f2D(:)-uu(:)))]
return
figure;

plot(xx,fval,'r-o',xx,func(xx),'b-<')
figure;
plot(xx,fval-func(xx),'r-o')
[norm(fval-func(xx)) max(abs(fval-func(xx)))]

% for 2D sparse grids solution

