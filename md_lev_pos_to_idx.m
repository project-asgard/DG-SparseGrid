function idx = md_lev_pos_to_idx(lev_vec,pos_vec,max_lev)

%%
% Given the lev and cell d-dim vectors, return the global element idx.

num_dimensions = numel(lev_vec);

assert(numel(pos_vec)==num_dimensions);

idx   = uint64(1);
stride = uint64(1);
for d=1:num_dimensions
    
    %%
    % ensure we are within the addressable space
    
    assert(lev_vec(d) <= max_lev);

    idx_1D  = lev_cell_to_1D_index(lev_vec(d),pos_vec(d));
    idx    = idx + (idx_1D-1)*stride;
    stride  = stride * 2^max_lev;
    
end

assert(idx >= 0);
assert(idx <= (uint64(2)^max_lev)^num_dimensions);

end