function fval = exact_solution_vector(HASHInv,pde,time)

% Returns the wavelet transformed exact solution

nDims = numel(pde.dimensions);

%%
% Loop over the number of dimensions (+time) to construct the analytic solution.

for d=1:nDims
    fList{d} = forwardMWT(pde.dimensions{d}.lev,pde.dimensions{d}.deg,...
        pde.dimensions{d}.domainMin,pde.dimensions{d}.domainMax,...
        pde.analytic_solutions_1D{d},pde.params,time);
end

ft = pde.analytic_solutions_1D{nDims+1}(time);
fval = combine_dimensions_D(fList,ft,HASHInv,pde);

end