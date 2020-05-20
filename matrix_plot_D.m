function [Meval,nodes] = matrix_plot_D(pde,dimension)

%%
% Generate the evaluation matrix and plotting points
% TODO : probably should be renamed to something more descriptive

%%
% Setup a few shortcuts

lev = dimension.lev;
deg = pde.deg;
xMin = dimension.domainMin;
xMax = dimension.domainMax;
FMWT = dimension.FMWT;

%%
% Jacobi of variable

n = 2^(lev);
h = (xMax-xMin)/n;
dof_1D = deg*n;
nodes = zeros(dof_1D,1);

%%
% Quadrature points (quad_x) and weights (quad_w)
[quad_x_interior_element,quad_w]=lgwt(deg,-1,1);

%%
% Add end points to plot

quad_x_left_element  = [-1 quad_x_interior_element']';
quad_x_right_element = [quad_x_interior_element' +1]';

dof = numel(quad_x_interior_element);

p_val       = lin_legendre(quad_x_interior_element,deg)*sqrt(1/h); % TODO : this call and normalization happens in several places. We should consolidate.
p_val_left  = lin_legendre(quad_x_left_element,deg)    *sqrt(1/h);
p_val_right = lin_legendre(quad_x_right_element,deg)   *sqrt(1/h);

Meval = sparse(dof*(n-2)+(dof+1)*2,dof_1D);

for i=0:n-1
    
    if i==0 
        quad_x = quad_x_left_element;
    elseif i==n-1
        quad_x = quad_x_right_element;
    else
        quad_x = quad_x_interior_element;
    end
    
    %%
    % Mapping from [-1,1] to physical domain
    xi = ((quad_x+1)/2+i)*h+xMin;
    
    %%
    % Coefficients for DG bases
    Iv = [deg*i+1:deg*(i+1)];
    if i==0
        Iu = [1:dof+1];
        Meval(Iu,Iv) = p_val_left;       
    elseif i==n-1
        Iu = [dof*i+1+1:dof*(i+1)+2];
        Meval(Iu,Iv) = p_val_right;       
    else
        Iu = [dof*i+1:dof*(i+1)]+1;
        Meval(Iu,Iv) = p_val;       
    end
    
    nodes(Iu) = xi;
    
end

%%
% Transform back to real space from wavelet space

Meval = Meval*FMWT';


    



