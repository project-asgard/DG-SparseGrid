function [child_elements_idx, cnt] = get_child_elements_idx (hash_table, idx, max_lev, refinement_method)

lev_vec = hash_table.elements.lev_p1(idx,:)-1;
pos_vec = hash_table.elements.pos_p1(idx,:)-1;

[child_elements_lev_vec, child_elements_pos_vec, cnt] = ...
    get_child_elements (lev_vec, pos_vec, max_lev, refinement_method);

child_elements_idx = [];
for i=1:cnt
    child_elements_idx(i) = lev_cell_to_element_index(...
        child_elements_lev_vec(i), child_elements_pos_vec(i), max_lev);
end

end