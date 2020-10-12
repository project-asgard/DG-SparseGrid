function moment_value = moment_integral(lev_vec, deg, coords_nD, fval_realspace, gfunc_nD, dimensions)

xmin = dimensions{1,1}.domainMin;
xmax = dimensions{1,1}.domainMax;
num_dim = length(dimensions);

Lev = lev_vec(1);
h = (xmax - xmin)/2^(Lev(1));

[quad_xx, quad_ww] = lgwt(deg, -1, 1);

quad_ww = 2^(-Lev)/2*quad_ww;

ww = repmat(quad_ww, 2^Lev(1), 1); 

if num_dim >= 2
    for i = 2:num_dim
         domainMin = dimensions{1,i}.domainMin;
         domainMax = dimensions{1,i}.domainMax;
         ww = kron(ww,ww)*(domainMax - domainMin);
    end
end

ww = ww.*(xmax - xmin);

[x, w] = lgwt(deg, 0, h);

this_dim_coord = coords_nD{1}(:);
jac = this_dim_coord.*0+1;
moment = this_dim_coord.*0+1;

% if num_dim == 1
%     moment = 
%     moment_value = sum(ww.*fval_realspace.*gfunc_nD(points).*dimensions{1}.jacobian(points));
% else
    for d=1:num_dim
        this_dim_coord = coords_nD{d}(:); 
        jac = jac .* dimensions{d}.jacobian(this_dim_coord);
        moment = moment .* gfunc_nD{d}(this_dim_coord);
    end
    moment_value = sum(ww.*fval_realspace.*moment.*jac);
% end
end
