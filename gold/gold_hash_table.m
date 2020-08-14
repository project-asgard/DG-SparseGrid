% Generate testing data for C++

addpath(genpath(pwd));



% element_table testing files
data_dir = strcat("generated-inputs", "/", "element_table", "/");
root = get_root_folder();
[stat,msg] = mkdir ([root,'/gold/',char(data_dir)]);

out_format = strcat(data_dir, "element_table_1_2_SG.dat");
num_dimensions = 1;
lev_vec = zeros(num_dimensions,1) + 2;
grid_type = 'SG';
[fwd1, inv1] = hash_table_nD(lev_vec, grid_type);
inv1_mat = cell2mat(inv1);
inv1_mat = reshape(inv1_mat,[3*num_dimensions,size(inv1,2)])';
coord = inv1_mat(:,1:2*num_dimensions);
filename = out_format;
write_octave_like_output(filename,coord);

out_format = strcat(data_dir, "element_table_2_3_SG.dat");
num_dimensions = 2;
lev_vec = zeros(num_dimensions,1) + 3;
grid_type = 'SG';
[fwd2, inv2] = hash_table_nD(lev_vec, grid_type);
inv2_mat = cell2mat(inv2);
inv2_mat = reshape(inv2_mat,[3*num_dimensions,size(inv2,2)])';
coord = inv2_mat(:,1:2*num_dimensions);
filename = out_format;
write_octave_like_output(filename,coord);

out_format = strcat(data_dir, "element_table_3_4_FG.dat");
num_dimensions = 3;
lev_vec = zeros(num_dimensions,1) + 4;
grid_type = 'FG';
[fwd3, inv3] = hash_table_nD(lev_vec, grid_type);
inv3_mat = cell2mat(inv3);
inv3_mat = reshape(inv3_mat,[3*num_dimensions,size(inv3,2)])';
coord = inv3_mat(:,1:2*num_dimensions);
filename = out_format;
write_octave_like_output(filename,coord);

%FIXME delete above
% element table tests

levels{1} = [9];
levels{2} = [5, 2];
levels{3} = [3, 2, 3];

out_format = strcat(data_dir, "table_%dd_%s.dat");
id_out_format = strcat(data_dir, "ids_%dd_%s.dat");

for i=1:size(levels, 2)
    
    [sparse_elements, sparse_elementsIDX] = hash_table_sparse_nD (levels{i}, max(levels{i}), 'SG');
    [level_row,level_col,level_val] = find(sparse_elements.lev_p1);
    [cell_row, cell_col, cell_val] = find(sparse_elements.pos_p1);
    num_entries = size(cell_val,1) / size(levels{i}, 2); 
    table = zeros(num_entries, size(levels{i},2) * 2);
    for j=1:num_entries
       id = sparse_elementsIDX(j);
       table(j, :) = horzcat(full(sparse_elements.lev_p1(id,:)), full(sparse_elements.pos_p1(id,:)));
    end
    
    
    [full_elements, full_elementsIDX] = hash_table_sparse_nD (levels{i}, max(levels{i}), 'FG');
    [cell_row, cell_col, cell_val] = find(full_elements.pos_p1);
    num_entries = size(cell_val,1) / size(levels{i}, 2); 
    full_table = zeros(num_entries, size(levels{i},2) * 2);
    for j=1:num_entries
        id = full_elementsIDX(j);
        full_table(j, :) = horzcat(full(full_elements.lev_p1(id,:)), full(full_elements.pos_p1(id,:)));
    end
        
    filename = sprintf(out_format, size(levels{i}, 2), "FG");
    write_octave_like_output(filename,full_table);
        
    filename = sprintf(out_format, size(levels{i}, 2), "SG");
    write_octave_like_output(filename,table);
    
    filename = sprintf(id_out_format, size(levels{i}, 2), "FG");
    write_octave_like_output(filename,full_elementsIDX');
    
    filename = sprintf(id_out_format, size(levels{i}, 2), "SG");
    write_octave_like_output(filename,sparse_elementsIDX');
    
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



