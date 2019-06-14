function [parentIdx] = get_my_parent_idx(pde,idx)

nDims = numel(pde.dimensions);

thisElemLevVec = pde.elements.lev_p1(idx,:)-1;
thisElemPosVec = pde.elements.pos_p1(idx,:)-1;

for d=1:nDims
    parentElemLevVec(d) = thisElemLevVec(d)-1;
    parentElemPosVec(d) = floor(thisElemPosVec(d)/2);  
end

parentIdx = lev_cell_to_element_index(pde,parentElemLevVec,parentElemPosVec);

end