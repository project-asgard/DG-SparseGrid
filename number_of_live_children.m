function [num_live_children, has_complete_children] = ...
    number_of_live_children (hash_table, lev_vec, pos_vec, max_lev, refinement_method)

assert (min(lev_vec)>=0);
assert (min(pos_vec)>=0);

%%
% get list of children

[child_elements_lev_vec, child_elements_pos_vec, num_children] = ...
    get_child_elements(lev_vec, pos_vec, max_lev, refinement_method);

%%
% check if children are live in the global list of elements

num_live_children = 0;
for ii=1:num_children
    this_element_lev_vec = child_elements_lev_vec(ii,:);
    this_element_pos_vec = child_elements_pos_vec(ii,:);
    this_idx = lev_cell_to_element_index(this_element_lev_vec,this_element_pos_vec,max_lev);
    if(hash_table.elements.type(this_idx) > 0)
        num_live_children = num_live_children + 1;
    end
end

%%
% mark this element if it is complete or incomplete

has_complete_children = false;
if num_live_children == num_children
    has_complete_children = true;
end
if num_children == 0
    has_complete_children = false;
end


end