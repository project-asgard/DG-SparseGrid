%function [Meval,node]=matrix_plot(lev,deg,xMin,xMax,FMWT)
function [Meval,nodes]=matrix_plot(dimension)

%%
% Generate the evaluation matrix and plotting points
% TODO : probably should be renamed to something more descriptive

%%
% Setup a few shortcuts

lev = dimension.lev;
deg = dimension.deg;
xMin = dimension.domainMin;
xMax = dimension.domainMax;
FMWT = dimension.FMWT;

%%
% Jacobi of variable

n = 2^(lev);
h = (xMax-xMin)/n;
dof_1D = deg*n;
nodes = zeros(dof_1D,1);
Meval = sparse(dof_1D,dof_1D);

%%
% Quadrature points (quad_x) and weights (quad_w)
[quad_x,quad_w]=lgwt(deg,-1,1);
p_val = legendre(quad_x,deg)*sqrt(1/h); % TODO : this call and normalization happens in several places. We should consolidate.

for i=0:n-1
    
    %%
    % Mapping from [-1,1] to physical domain
    xi = ((quad_x+1)/2+i)*h+xMin;
    
    %%
    % Coefficients for DG bases
    Iu = [deg*i+1:deg*(i+1)];
    
    Meval(Iu,Iu) = p_val;
    nodes(Iu) = xi;
    
end

%%
% Transform back to real space from wavelet space

Meval = Meval*FMWT';


    



