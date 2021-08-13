function [mass] = calculate_mass(pde,opts,coord,fval_realspace)

num_dims = numel(pde.dimensions);

mass_func = @(x,p,t) x.*0+1;
for d=1:num_dims
    moment_func_nD{d} = mass_func;
end

mass = moment_integral(pde.get_lev_vec,opts.deg,coord,fval_realspace,moment_func_nD,pde.dimensions);

end