function fval = initial_condition_vector(pde,opts,hash_table,t)

% Returns the multi-D wavelet space initial condition

if opts.many_solution_capable
    
    fval = md_eval_function(opts,opts.deg,pde.dimensions,pde.params, ...
        pde.initial_conditions,hash_table,pde.transform_blocks,t);
    
else
    
    dims = pde.dimensions;
    nDims = numel(dims);
    
    for d=1:nDims
        fList{d} = forward_wavelet_transform(opts.deg,dims{d}.lev,...
            dims{d}.min,dims{d}.max,...
            dims{d}.init_cond_fn,pde.params,pde.transform_blocks, t);
    end
    
    ft = 1;
    fval = combine_dimensions_D(opts.deg, fList, ft, ...
        hash_table, opts.use_oldhash);
    
end

end
