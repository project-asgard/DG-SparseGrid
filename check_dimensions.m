function pde = check_dimensions(pde,opts)

dims = pde.dimensions;

nDims = numel(dims);

for d=1:nDims
   
    pde.dimensions{d} = check_dimension(opts,nDims,dims{d});
    
end

end
