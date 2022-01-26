function fval = initial_condition_vector(pde,opts,hash_table,t)

% Returns the multi-D wavelet space initial condition

num_sol = numel(pde.initial_conditions);

fval = cell(num_sol,1);

for i = 1 : num_sol

    fval{i} = md_eval_function(opts,opts.deg,pde.dimensions,pde.params, ...
      pde.initial_conditions{i},hash_table,pde.transform_blocks,t);
  
end

end
