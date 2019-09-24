function [parentIdx] = get_my_parent_idx(num_dims, hash_table, idx, max_lev)

thisElemLevVec = hash_table.elements.lev_p1(idx,:)-1;
thisElemPosVec = hash_table.elements.pos_p1(idx,:)-1;

for d=1:num_dims
    parentElemLevVec(d) = thisElemLevVec(d)-1;
    parentElemPosVec(d) = floor(thisElemPosVec(d)/2);  
end

parentIdx = lev_cell_to_element_index(parentElemLevVec,parentElemPosVec, max_lev);

end