function [mass] = calculate_mass(pde,opts,fval,hash_table)

num_dims = numel(pde.dimensions);

mass_func = @(x,p,t) x.*0+1;
for d=1:num_dims+1
    moment_func_nD{d} = mass_func;
end

mass_bool = 0; %Want the mass matrix applied
M_one = md_eval_function(opts,opts.deg,pde.dimensions,pde.params,...
    {moment_func_nD},hash_table, pde.transform_blocks,0,mass_bool);

mass = M_one'*fval;
end