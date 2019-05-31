function fval = initial_condition_vector(HASHInv,pde,time)

dims = pde.dimensions;
nDims = numel(dims);

for d=1:nDims
    fList{d} = forwardMWT(pde,d,dims{d}.init_cond_fn,time);
end

ft = 1;
fval = combine_dimensions_D(fList,ft,HASHInv,pde);

end