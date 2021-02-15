function ans = moment_integral(lev_vec,deg,f,md_gfunc,dims,nodes,dims_subset_idx)

% evaluation of some function (md_gfunc) moment of a real-space distribution function f

%assert(numel(squeeze(size(f))) == numel(dims)); % this function now accepts the m x n (for 2D) shaped f

num_dims_from_f = numel(size(f));
if num_dims_from_f <= 2
   if numel(f(:)) == max(size(f))
       num_dims_from_f = 1;
   end
end

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
    dim_min = dims{1,d1}.min;
    dim_max = dims{1,d1}.max;
    if ~isempty(find(vecdim==d1)) % only apply quadrature weights for those dimensions we are integrating over
        ww = kron(ww,ww_1D)*(dim_max - dim_min);
    else
        ww = kron(ww,ww_1D.*0+1);
    end
end

ww = reshape(ww,size(f));

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


% we want the following (which can be summed over to give the integral due
% to the use the quadrature weights ww at the quadrature points xx, yy 
% (2D example)
% ww(xx,yy) * f(xx,yy) * g(xx,yy) * J(xx,yy)
% which can be written as the product of 1D funcions as
% 1 * ww(xx) * f(xx) * g(xx) * J(xx) * ww(yy) * f(yy) * g(yy) * J(yy)
% where xx and yy are of md_ type. 

% e.g., 2D, r, th (spherical with ph integrated away)

% dl = gr(r)gr(th) r_hat + gth(r)gth(th) th_hat
% gr(r) = 1, gr(th) = 1
% gth(r) = r, gth(th) = 1

if num_dims > 1    
    f = permute(f,flip(1:num_dims)); % again, because f dims are in the wrong order at the main level :(
end

for d1 = 1:num_dims
    this_dim_coord = coords{d1};
    if ~isempty(find(vecdim==d1)) % only apply moment func and jac for those dimensions we are integrating over
        moment = moment .* md_gfunc{d1}(this_dim_coord);
        
        for d2 = 1:num_dims
            this_dim_coord = coords{d2};
            jac = jac .* dims{d1}.jacobian(this_dim_coord);
        end
    end
end
moment = reshape(moment,size(f));
jac = reshape(moment, size(f));
md_moment = ww.*f.*moment.*jac;
ans = sum(md_moment,vecdim);

end
