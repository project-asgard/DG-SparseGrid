function fval = exact_solution_vector(pde,opts,hash_table,time)

% Returns the wavelet transformed exact solution

nDims = numel(pde.dimensions);

%%
% Loop over the number of dimensions (+time) to construct the analytic solution.

for d=1:nDims
    fList{d} = forward_wavelet_transform(pde,d,pde.analytic_solutions_1D{d},time);
end

ft = pde.analytic_solutions_1D{nDims+1}(time);
fval = combine_dimensions_D(pde,opts,fList,ft,hash_table);

end