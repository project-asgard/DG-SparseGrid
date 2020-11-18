function fval = source_vector(pde,opts,hash_table,t)

% Returns the wavelet transformed source

fval = md_eval_function(opts, opts.deg, pde.dimensions, ...
    pde.params, pde.sources, hash_table, pde.transform_blocks, t);

end
