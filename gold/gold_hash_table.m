% Generate testing data for C++

addpath(genpath(pwd));

% element_table testing files
data_dir = strcat("generated-inputs", "/", "element_table", "/");
root = get_root_folder();
[stat,msg] = mkdir ([root,'/gold/',char(data_dir)]);

% element table tests
levels = {};
levels{1} = [7];
levels{2} = [5, 2];
levels{3} = [3, 2, 3];
max_adapt_lev = 8;
out_format = strcat(data_dir, "table_%dd_%s.dat");
id_out_format = strcat(data_dir, "ids_%dd_%s.dat");
child_id_out_format = strcat(data_dir, "child_ids_%dd_%s.dat");
for i=1:size(levels, 2)
    [sparse_elements, sparse_elementsIDX] = hash_table_sparse_nD (levels{i}, max_adapt_lev, 'SG');
    [level_row,level_col,level_val] = find(sparse_elements.lev_p1);
    [cell_row, cell_col, cell_val] = find(sparse_elements.pos_p1);
    num_entries = size(cell_val,1) / size(levels{i}, 2); 
    table = zeros(num_entries, size(levels{i},2) * 2);
    for j=1:num_entries
       id = sparse_elementsIDX(j);
       table(j, :) = horzcat(full(sparse_elements.lev_p1(id,:)), full(sparse_elements.pos_p1(id,:))) - 1;
    end
    
    
    [full_elements, full_elementsIDX] = hash_table_sparse_nD (levels{i}, max_adapt_lev, 'FG');
    [cell_row, cell_col, cell_val] = find(full_elements.pos_p1);
    num_entries = size(cell_val,1) / size(levels{i}, 2); 
    full_table = zeros(num_entries, size(levels{i},2) * 2);
    for j=1:num_entries
        id = full_elementsIDX(j);
        full_table(j, :) = horzcat(full(full_elements.lev_p1(id,:)), full(full_elements.pos_p1(id,:))) - 1;
    end
        
    filename = sprintf(out_format, size(levels{i}, 2), "FG");
    write_octave_like_output(filename,full_table);
        
    filename = sprintf(out_format, size(levels{i}, 2), "SG");
    write_octave_like_output(filename,table);
    
    filename = sprintf(id_out_format, size(levels{i}, 2), "FG");
    write_octave_like_output(filename,(full_elementsIDX'-1));
    
    filename = sprintf(id_out_format, size(levels{i}, 2), "SG");
    write_octave_like_output(filename,(sparse_elementsIDX'-1));
    
    
    sparse_children = generate_child_ids(sparse_elementsIDX, size(levels{i}, 2), ...
                                         max_adapt_lev);
    full_children = generate_child_ids(full_elementsIDX, size(levels{i}, 2), ...
                                         max_adapt_lev);
                                     
    filename = sprintf(child_id_out_format, size(levels{i}, 2), "SG");
    write_octave_like_output(filename,(sparse_children-1));

    filename = sprintf(child_id_out_format, size(levels{i}, 2), "FG");
    write_octave_like_output(filename,(full_children-1));
    
end

% 1d indexing tests
pairs{1} = [0, 0];
pairs{2} = [1, 0];
pairs{3} = [2, 1];
pairs{4} = [12, 5];
pairs{5} = [7, 0];
pairs{6} = [4, 6];
pairs{7} = [9, 3];
pairs{8} = [30, 20];
out_format = strcat(data_dir, "1d_index_%d_%d.dat");
for i=1:size(pairs, 2)
    id = lev_cell_to_1D_index(pairs{i}(1), pairs{i}(2));
    filename = sprintf(out_format, pairs{i}(1), pairs{i}(2));
    write_octave_like_output(filename, id);
end


function [child_ids] = generate_child_ids(parent_element_ids, num_dims, max_lev)
    
    child_ids = [];
    refinement_method = 1; % this is the only method implemented in C++ currently
    for i=1:length(parent_element_ids)
        children_i = get_child_elements_idx (num_dims, max_lev,...
                    parent_element_ids(i), refinement_method);
        child_ids = [child_ids ; children_i];
    end
    
    
end
