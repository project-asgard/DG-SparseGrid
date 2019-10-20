function fval = initial_condition_vector(pde, opts, hash_table, t)

% Returns the multi-D wavelet space initial condition

if opts.many_solution_capable
    
    fval = md_eval_function(opts, pde.deg, pde.dimensions, pde.params, ...
        pde.initial_conditions, hash_table, t);
    
else
    
    dims = pde.dimensions;
    nDims = numel(dims);
    
    for d=1:nDims
        fList{d} = forward_wavelet_transform(pde.deg,dims{d}.lev,...
            dims{d}.domainMin,dims{d}.domainMax,...
            dims{d}.init_cond_fn,pde.params,t);
    end
    
    ft = 1;
    fval = combine_dimensions_D(pde.deg, fList, ft, ...
        hash_table, opts.use_oldhash);
    
end

end