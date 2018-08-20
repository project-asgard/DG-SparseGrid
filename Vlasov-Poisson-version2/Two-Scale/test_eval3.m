% Test
% Line 3-36 is to set up fval, which can be removed for the real code
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

% define two functions and func = func_x*func_y
func_x = @(x,params)(sin(pi*x));
func_y = @(y,params)(sin(pi*y));


% compute the 1D coefficient 
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
% Method :: evaluate every basis at the given points (xx) and (vv)
% fx_loc = phi_x (xx) with size DoF1DxNum_xx
% fy_loc = phi_y (yy) with size DoF1DxNum_vv
% Loop over every point with (xx_i,yy_j) and combine the value
%   for fx_loc(:,i) and fy_loc(:,j)

% if we have computed fval, we can just start with following
% Input:: Lstart, Lend, Lev, Deg, xx, yy
xx = [Lstart:0.1:Lend]';
yy = [Lstart:0.1:Lend]';

[fx_loc] = EvalWavPoint4(Lstart,Lend,Lev,Deg,xx);
[fy_loc] = EvalWavPoint4(Lstart,Lend,Lev,Deg,yy);


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


% Plot for varification
[XX,YY] = meshgrid(xx,yy);
nx = length(xx);ny = length(yy);
uu = func_x(XX,params).*func_y(YY,params);


figure(1);
subplot(1,3,1)
surf(XX,YY,f2D);
title('Solution on Given Points')
subplot(1,3,2)
surf(XX,YY,uu);
title('Exact Solution')


[norm(f2D(:)-uu(:)) max(abs(f2D(:)-uu(:)))]


% for new points, we shall use interpolation to evaluate
xnew = [Lstart:0.05:Lend]';
ynew = [Lstart:0.05:Lend]';

[XXnew,YYnew] = meshgrid(xnew,ynew);
nxnew = length(xnew);
nynew = length(ynew);

f2D_inter = interp2(XX,YY,full(f2D),XXnew,YYnew);
figure(1);
subplot(1,3,3)
surf(XXnew,YYnew,f2D_inter);
title('Interpolation')
