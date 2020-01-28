function fval = exact_solution_vector(pde,opts,hash_table,t)

% Returns the wavelet transformed exact solution

if opts.many_solution_capable
    
    fval = md_eval_function(opts, pde.deg, pde.dimensions, ...
        pde.params, pde.solutions, hash_table, t);
    
else
    
    nDims = numel(pde.dimensions);
    
    %%
    % Loop over the number of dimensions (+time) to construct the analytic solution.
    
    for d=1:nDims
        fList{d} = forward_wavelet_transform(pde.deg,pde.dimensions{d}.lev,...
            pde.dimensions{d}.domainMin,pde.dimensions{d}.domainMax,...
            pde.analytic_solutions_1D{d},pde.params,t);
    end
    
    ft = pde.analytic_solutions_1D{nDims+1}(t);
    fval = combine_dimensions_D(pde.deg,fList,ft,hash_table,opts.use_oldhash);
    
end

end