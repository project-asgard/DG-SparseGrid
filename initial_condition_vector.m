function fval = initial_condition_vector(pde,opts,hash_table,t)

% Returns the multi-D wavelet space initial condition

fval = md_eval_function(opts,opts.deg,pde.dimensions,pde.params, ...
    pde.initial_conditions,hash_table,pde.transform_blocks,t);

end
