function fval = exact_solution_vector(pde,opts,hash_table,t)

% Returns the wavelet transformed exact solution

fval = md_eval_function(opts, opts.deg, pde.dimensions, ...
    pde.params, pde.solutions, hash_table, pde.transform_blocks, t);

end
