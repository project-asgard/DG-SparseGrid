function pde = checkAllDimensions(pde)

dims = pde.dimensions;

nDims = numel(dims);

for d=1:nDims
   
    pde.dimensions{d} = checkDimension(nDims,dims{d});
    
end

end