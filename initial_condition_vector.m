function fval = initial_condition_vector(pde, opts, hash_table, time)

dims = pde.dimensions;
nDims = numel(dims);

for d=1:nDims
    fList{d} = forwardMWT(pde,d,dims{d}.init_cond_fn,time);
end

ft = 1;
fval = combine_dimensions_D(pde, opts, fList, ft, hash_table);

end