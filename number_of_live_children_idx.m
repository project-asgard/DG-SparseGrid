function [num_live_children, has_complete_children] = ...
    number_of_live_children_idx (hash_table, idx, max_lev, refinement_method)

lev_vec = hash_table.elements.lev_p1(idx,:);
pos_vec = hash_table.elements.pos_p1(idx,:);


[num_live_children, has_complete_children] = ...
    number_of_live_children (hash_table, lev_vec, pos_vec, max_lev, refinement_method);

end