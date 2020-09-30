function [Meval,nodes] = matrix_plot_D(pde,opts,dimension)

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

if opts.uniform_output
    % output on uniformly spaced grid
    % note that this uses the left of the dual value options
    N = deg+2;
    delta = 2/N;
    quad_x_interior_element = linspace(-1,1-delta,N)';
    quad_x_left_element  = quad_x_interior_element;
    quad_x_right_element = [quad_x_interior_element' +1]';
    dof = numel(quad_x_interior_element);
    Meval = sparse(dof*(N-1)+(dof+1)*1,dof_1D);
elseif opts.add_end_points
    % output on quadrature points (quad_x) with end points
    [quad_x_interior_element,~]=lgwt(deg,-1,1);
    quad_x_left_element  = [-1 quad_x_interior_element']';
    quad_x_right_element = [quad_x_interior_element' +1]';
    dof = numel(quad_x_interior_element);
    Meval = sparse(dof*(n-2)+(dof+1)*2,dof_1D);
else
    % output on quadrature points (quad_x)
    [quad_x_interior_element,~]=lgwt(deg,-1,1);
    dof = numel(quad_x_interior_element);
    Meval = sparse(dof_1D,dof_1D);
end

p_val       = lin_legendre(quad_x_interior_element,deg)*sqrt(1/h); % TODO : this call and normalization happens in several places. We should consolidate.
if opts.add_end_points || opts.uniform_output
    p_val_left  = lin_legendre(quad_x_left_element,deg)    *sqrt(1/h);
    p_val_right = lin_legendre(quad_x_right_element,deg)   *sqrt(1/h);
end

for i=0:n-1
    
    quad_x = quad_x_interior_element;
    
    if opts.add_end_points || opts.uniform_output    
        if i==0
            quad_x = quad_x_left_element;
        elseif i==n-1
            quad_x = quad_x_right_element;
        end
    end
    
    %%
    % Mapping from [-1,1] to physical domain
    xi = ((quad_x+1)/2+i)*h+xMin;
    
    %%
    % Coefficients for DG bases
    Iv = [deg*i+1:deg*(i+1)];
    if opts.uniform_output
        if i==0
            Iu = [1:dof];
            Meval(Iu,Iv) = p_val_left;
        elseif i==n-1
            Iu = [dof*i+1:dof*(i+1)+1];
            Meval(Iu,Iv) = p_val_right;
        else
            Iu = [dof*i+1:dof*(i+1)];
            Meval(Iu,Iv) = p_val;
        end
    elseif opts.add_end_points
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
    else      
        Iu = [dof*i+1:dof*(i+1)];
        Meval(Iu,Iv) = p_val;
    end
    
    nodes(Iu) = xi;
    
end

%%
% Transform back to real space from wavelet space

Meval = Meval*FMWT';

