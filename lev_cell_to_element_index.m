function eIdx = lev_cell_to_element_index(pde,lev,pos)

%%
% Given the lev and cell d-dim vectors, return the global element idx.

dims  = pde.dimensions;
nDims = numel(dims);

assert(numel(lev)==nDims);
assert(numel(pos)==nDims);

idx1D = lev_cell_to_singleD_index(lev,pos);

eIdx   = uint64(0);
stride = uint64(1);
for d=1:nDims
    dim = dims{d};
    eIdx = eIdx + idx1D(d)*stride;
    assert(dim.lev <= pde.maxLev);
    stride = stride * 2^pde.maxLev;
end

assert(eIdx >= 0);

end