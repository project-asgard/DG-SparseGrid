function pde = check_dimensions(pde)

dims = pde.dimensions;

nDims = numel(dims);

for d=1:nDims
   
    pde.dimensions{d} = check_dimension(nDims,dims{d});
    pde.lev_vec(d) = pde.dimensions{d}.lev;
    
end

end