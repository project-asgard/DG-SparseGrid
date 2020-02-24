function [parent_elements_idx, num_parents] = get_parent_elements_idx (hash_table, idx, max_lev, refinement_method)

lev_vec = hash_table.elements.lev_p1(idx,:)-1;
pos_vec = hash_table.elements.pos_p1(idx,:)-1;

[parent_elements_lev_vec, parent_elements_pos_vec, num_parents] = ...
    get_parent_elements (lev_vec, pos_vec, max_lev, refinement_method);

if num_parents > 0
    for i=1:num_parents
        parent_elements_idx(i) = lev_cell_to_element_index(...
            parent_elements_lev_vec(i,:), parent_elements_pos_vec(i,:), max_lev);
    end
else
    parent_elements_idx = NaN;
end

end