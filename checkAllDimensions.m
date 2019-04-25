function pde = checkAllDimensions(pde)

dims = pde.dimensions;

nDims = numel(dims);

for d=1:nDims
   
    pde.dims{d} = checkDimension(nDims,dims{d});
    
end

end