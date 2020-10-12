function [Meval,nodes,nodes_count] = matrix_plot_D(pde,opts,dimension,input_nodes)

%%
% Generate the evaluation matrix and plotting points
% TODO : probably should be renamed to something more descriptive

user_provided_nodes = false;
if nargin>3
    user_provided_nodes = true;
end

%%
% Setup a few shortcuts

lev = dimension.lev;
deg = opts.deg;
xMin = dimension.domainMin;
xMax = dimension.domainMax;
% FMWT = dimension.FMWT;

%%
% Jacobi of variable

n = 2^(lev);
h = (xMax-xMin)/n;
dof_1D = deg*n;

if user_provided_nodes
    % output at a user specified set of locations
    input_nodes = reshape(input_nodes,numel(input_nodes),1);
    num_nodes = numel(input_nodes);
    nodes_count = input_nodes .* 0;
    dof = num_nodes*deg;
    Meval = sparse(num_nodes,dof_1D);
elseif strcmp(opts.output_grid,'uniform')
    % output on uniformly spaced grid
    % note that this uses the left of the dual value options
    N = deg+2;
    delta = 2/N;
    quad_x_interior_element = linspace(-1,1-delta,N)';
    quad_x_left_element  = quad_x_interior_element;
    quad_x_right_element = [quad_x_interior_element' +1]';
    dof = numel(quad_x_interior_element);
    Meval = sparse(dof*(N-1)+(dof+1)*1,dof_1D);
elseif strcmp(opts.output_grid,'quadrature_with_end_points')
    % output on quadrature points (quad_x) with end points
    % note that this will produce dual valued output
    [quad_x_interior_element,~]=lgwt(deg,-1,1);
    quad_x_left_element  = [-1 quad_x_interior_element']';
    quad_x_right_element = [quad_x_interior_element' +1]';
    dof = numel(quad_x_interior_element);
    Meval = sparse(dof*(n-2)+(dof+1)*2,dof_1D);
else
    % output on quadrature points (quad_x) without end points
    [quad_x_interior_element,~]=lgwt(deg,-1,1);
    dof = numel(quad_x_interior_element);
    Meval = sparse(dof_1D,dof_1D);
end

cnt = 0;
for i=0:n-1
    
    if user_provided_nodes
        % which nodes are in this element
        quad_x = ((input_nodes - xMin)/h - i)*2-1;
        in_this_node = quad_x>=-1 & quad_x<=+1;
        nodes_count(in_this_node) = nodes_count(in_this_node)+1; % for averaging dual values later
        quad_x = quad_x(in_this_node);
    else
        quad_x = quad_x_interior_element;
        if strcmp(opts.output_grid,'quadrature_with_end_points') || strcmp(opts.output_grid,'uniform')
            if i==0
                quad_x = quad_x_left_element;
            elseif i==n-1
                quad_x = quad_x_right_element;
            end
        end
    end
    
    if ~isempty(quad_x)
        
        p_val = lin_legendre(quad_x,deg)*sqrt(1/h); % TODO : this call and normalization happens in several places. We should consolidate.
        
        %%
        % Mapping from [-1,1] to physical domain
        xi = ((quad_x+1)/2+i)*h+xMin;
        
        if user_provided_nodes
            assert(norm(xi-input_nodes(in_this_node))<=1e-10);
        end
        
        %%
        % Coefficients for DG bases
        Iv = [deg*i+1:deg*(i+1)];
        if user_provided_nodes
            Iu = [cnt+1:cnt+numel(quad_x)];
            cnt = cnt + numel(quad_x);
        elseif strcmp(opts.output_grid,'uniform')
            if i==0
                Iu = [1:dof];
            elseif i==n-1
                Iu = [dof*i+1:dof*(i+1)+1];
            else
                Iu = [dof*i+1:dof*(i+1)];
            end
        elseif strcmp(opts.output_grid,'quadrature_with_end_points')
            if i==0
                Iu = [1:dof+1];
            elseif i==n-1
                Iu = [dof*i+1+1:dof*(i+1)+2];
            else
                Iu = [dof*i+1:dof*(i+1)]+1;
            end
        else
            Iu = [dof*i+1:dof*(i+1)];
        end
        
        Meval(Iu,Iv) = p_val;
        nodes(Iu) = xi;
        
    else
        error('ERROR : quad_x is empty');
    end
end

%%
% Transform back to real space from wavelet space

%Meval = Meval*FMWT';
trans_side = 'RT';
Meval = apply_FMWT_blocks(lev, pde.transform_blocks, Meval, trans_side); 

