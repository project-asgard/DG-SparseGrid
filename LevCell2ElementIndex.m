function eIdx = LevCell2ElementIndex(pde,lev,cell)

%%
% Given the lev and cell d-dim vectors, return the global element idx.

dims = pde.dimensions;
nDims = numel(dims);

assert(numel(lev)==nDims);
assert(numel(cell)==nDims);

idx1D = LevCell2index(lev,cell);

eIdx = 0;
stride = 1;
for d=1:nDims
    dim = dims{d};
    eIdx = eIdx + idx1D(d)*stride;
    stride = stride * 2^dim.lev;
end

end