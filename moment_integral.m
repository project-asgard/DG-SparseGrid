function ans = moment_integral(lev_vec,deg,f,md_gfunc,dims,nodes,dims_subset_idx)

% evaluation of some function (md_gfunc) moment of a real-space distribution function f

assert(numel(size(f)) == numel(dims)); % this function now accepts the m x n (for 2D) shaped f

num_dims = length(dims);

% select which dimensions to do the moment over
vecdim = 1:num_dims;
if nargin>6
    if ~isempty(dims_subset_idx)
        vecdim = dims_subset_idx;
    end
end

lev = lev_vec(1);
if numel(lev_vec)>1
assert(sum(lev_vec./lev)==num_dims); % for now assert that lev is the same for all dims
end

% this implementation assumes that the solution is provide on the
% original quadrature points, i.e., asgard('output_grid','quadrature'),
% which is the default, and there is a catch to assert this in OPTS.m

% get quadrature points and weights
[quad_xx, quad_ww] = lgwt(deg, -1, 1);

% scale the weights becuase why? Lin?
quad_ww = 2^(-lev)/2*quad_ww;

ww_1D = repmat(quad_ww, 2^lev(1), 1);
ww = 1;
for d1 = 1:num_dims
    min = dims{1,d1}.min;
    max = dims{1,d1}.max;
    if ~isempty(find(vecdim==d1)) % only apply quadrature weights for those dimensions we are integrating over
        ww = kron(ww,ww_1D)*(max - min);
    else
        ww = kron(ww,ww_1D.*0+1);
    end
end

ww = reshape(ww,size(f))';

jac = ones('like',f);
moment = ones('like',f);

% create the md grid (because coords at the main level is in the wrong order)
if num_dims == 1
    coords{1} = nodes{1};
elseif num_dims == 2
    [coords{1},coords{2}] = ndgrid(nodes{1},nodes{2});
elseif num_dims ==3
    [coords{1},coords{2},coords{3}] = ndgrid(nodes{1},nodes{2},nodes{3});
else
    error('num_dims>3 not supported');
end

this_dim_coord = coords_nD{1}(:);
jac = this_dim_coord.*0+1;
moment = this_dim_coord.*0+1;

for d=1:num_dim
    this_dim_coord = coords_nD{d}(:);
    jac = jac .* dimensions{d}.moment_dV(this_dim_coord);
    moment = moment .* gfunc_nD{d}(this_dim_coord);
end

end
