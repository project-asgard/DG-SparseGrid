function [child_elements_idx, cnt] = get_child_elements_idx (num_dims, max_lev, idx, refinement_method)

% lev_vec = hash_table.elements.lev_p1(idx,:)-1;
% pos_vec = hash_table.elements.pos_p1(idx,:)-1;

[lev_vec, pos_vec] = md_idx_to_lev_pos(num_dims, max_lev, idx);

[child_elements_lev_vec, child_elements_pos_vec, cnt] = ...
    get_child_elements (lev_vec, pos_vec, max_lev, refinement_method);

child_elements_idx = zeros(cnt,1);
for i=1:cnt
    child_elements_idx(i) = md_lev_pos_to_idx(...
        child_elements_lev_vec(i,:), child_elements_pos_vec(i,:), max_lev);
end

end