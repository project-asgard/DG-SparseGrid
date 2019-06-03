function eIdx = lev_cell_to_element_index(pde,lev,pos)

%%
% Given the lev and cell d-dim vectors, return the global element idx.

dims  = pde.dimensions;
num_dimensions = numel(dims);

assert(numel(lev)==num_dimensions);
assert(numel(pos)==num_dimensions);

eIdx   = uint64(1);
stride = uint64(1);
for d=1:num_dimensions
    
    
    %%
    % ensure we are within the addressable space
    
    dim = dims{d};
    assert(dim.lev <= pde.maxLev);

    idx_1D  = lev_cell_to_singleD_index(lev(d),pos(d));
%     fprintf('dim: %i, idx: %i\n',d,idx_1D);
    eIdx    = eIdx + (idx_1D-1)*stride;
    stride  = stride * 2^pde.maxLev;
    
end

assert(eIdx >= 0);
assert(eIdx <= (uint64(2)^pde.maxLev)^num_dimensions);

end