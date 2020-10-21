function moment_value = moment_integral(lev_vec, deg, coords_nD, fval_realspace, gfunc_nD, dimensions)

xmin = dimensions{1,1}.min;
xmax = dimensions{1,1}.max;
num_dim = length(dimensions);

Lev = lev_vec(1);
h = (xmax - xmin)/2^(Lev(1));

[quad_xx, quad_ww] = lgwt(deg, -1, 1);

quad_ww = 2^(-Lev)/2*quad_ww;

ww_1D = repmat(quad_ww, 2^Lev(1), 1); 
ww = ww_1D;

if num_dim >= 2
    for i = 2:num_dim
         min = dimensions{1,i}.min;
         max = dimensions{1,i}.max;
         ww = kron(ww,ww)*(max - min);
    end
end

ww = ww.*(xmax - xmin);

this_dim_coord = coords_nD{1}(:);
jac = this_dim_coord.*0+1;
moment = this_dim_coord.*0+1;

for d=1:num_dim
    this_dim_coord = coords_nD{d}(:);
    jac = jac .* dimensions{d}.jacobian(this_dim_coord);
    moment = moment .* gfunc_nD{d}(this_dim_coord);
end
moment_value = sum(ww.*fval_realspace.*moment.*jac);

end
