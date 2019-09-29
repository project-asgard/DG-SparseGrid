function [parent_elements_idx] = get_element_parents (lev_vec, pos_vec, max_lev, refinement_method)

num_dims = numel(lev_vec);

if refinement_method == 1
    
    for d=1:num_dims
        
        above_lev_vec = lev_vec;
        above_lev_vec(d) = above_lev_vec(d)-1;
        above_pos_vec = pos_vec;
        above_pos_vec(d) = floor(above_pos_vec(d)/2);
        above_idx = lev_cell_to_element_index(above_lev_vec, above_pos_vec, max_lev);
        
        parent_elements_idx(d) = above_idx;
        
    end
    
else
    
    error('refinement_method ~= 1 not implemented yet');
    
end

end