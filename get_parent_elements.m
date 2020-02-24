function [parent_elements_lev_vec, parent_elements_pos_vec, num_parents] = ...
    get_parent_elements (lev_vec, pos_vec, max_lev, refinement_method)

num_dims = numel(lev_vec);

if refinement_method == 1
    
    for d=1:num_dims
        
        parent_lev_vec = lev_vec;
        parent_lev_vec(d) = parent_lev_vec(d)-1;
        parent_pos_vec = pos_vec;
        parent_pos_vec(d) = floor(parent_pos_vec(d)/2);
        
        %         above_idx = lev_cell_to_element_index(above_lev_vec, above_pos_vec, max_lev);
        %
        %         parent_elements_idx(d) = above_idx;
        
        if min(parent_lev_vec) >=0 % do not add parents above lev=0
            
            parent_elements_lev_vec(d,:) = parent_lev_vec;
            parent_elements_pos_vec(d,:) = parent_pos_vec;
            
        end
    end
    
else
    
    error('refinement_method ~= 1 not implemented yet');
    
end

[num_parents,num_dims2] = size(parent_elements_lev_vec);
assert(num_dims == num_dims2);

end