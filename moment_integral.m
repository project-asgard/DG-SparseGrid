function moment_value = moment_integral(lev_vec, deg, coords_nD, fval_realspace, gfunc_nD, dimensions, subset_dimensions)

xmin = dimensions{1,1}.min;
xmax = dimensions{1,1}.max;
num_dim = length(dimensions);

Lev = lev_vec(1);
h = (xmax - xmin)/2^(Lev(1));
num_elements = 2^Lev*deg;
submag = length(fval_realspace{1,1}).^subset_dimensions;
if subset_dimensions > 0
    for i = 1:submag
        for j = 1:num_elements
            for k = 1:num_elements
                sub_fval_realspace(i,:) = fval_realspace{1,1}(j,k,:);
            end
        end
    end
    fval_realspace = sub_fval_realspace;
else
    mag = length(fval_realspace{1,1}).^num_dim;
    fval_realspace = reshape(fval_realspace{1,1}, mag, 1);
end


[quad_xx, quad_ww] = lgwt(deg, -1, 1);

quad_ww = 2^(-Lev)/2*quad_ww;

ww_1D = repmat(quad_ww, 2^Lev(1), 1); 
ww = ww_1D;

if subset_dimensions >= 2 
    for i = 2:subset_dimensions
         min = dimensions{1,i}.min;
         max = dimensions{1,i}.max;
         ww = kron(ww,ww_1D)*(max - min);
    end
elseif num_dim >=2
   for i = 2:num_dim
         min = dimensions{1,i}.min;
         max = dimensions{1,i}.max;
         ww = kron(ww,ww_1D)*(max - min);
    end
end

ww = ww.*(xmax - xmin);
if subset_dimensions > 0
    for i = 1:submag
        for j = 1:num_elements
            for k = 1:num_elements
                this_dim_coord(i,:) = coords_nD{1,1}(j,k,:);
            end
        end
    end
else
    this_dim_coord = coords_nD{1}(:);
end
    jac = this_dim_coord.*0+1;
    moment = this_dim_coord.*0+1;

this_dim_coord = coords_nD{1}(:);
jac = this_dim_coord.*0+1;
moment = this_dim_coord.*0+1;

for d=1:num_dim
    this_dim_coord = coords_nD{d}(:);
    jac = jac .* dimensions{d}.moment_dV(this_dim_coord);
    moment = moment .* gfunc_nD{d}(this_dim_coord);
end
      
moment_value = sum(ww.*fval_realspace.*moment.*jac);

end
