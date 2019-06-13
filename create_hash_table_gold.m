% Generate testing data for C++ 

addpath(genpath(pwd));

% element_table testing files
element_dir = strcat(pwd, "/", "generated-inputs", "/", "element_table", "/");
[status,msg] = mkdir (element_dir);

out_format = strcat(element_dir, "element_table_1_1_SG.dat");
num_dimensions = 1;
lev_vec = zeros(num_dimensions,1) + 1;
grid_type = 'SG';
[fwd1, inv1] = create_hash_table(lev_vec, grid_type);
% for i=1:size(inv1,2)
    inv1_mat = cell2mat(inv1);
    inv1_mat = reshape(inv1_mat,[3*num_dimensions,size(inv1,2)])';
    coord = inv1_mat(:,1:2*num_dimensions);
%   coord = inv1{i}(1:2*num_dimensions);
%   filename = sprintf(out_format, i);
    filename = out_format;
  write_octave_like_output(filename,coord);
% end

out_format = strcat(element_dir, "element_table_2_3_SG.dat");
num_dimensions = 2;
lev_vec = zeros(num_dimensions,1) + 3;
grid_type = 'SG';
[fwd2, inv2] = create_hash_table(lev_vec, grid_type);
    inv2_mat = cell2mat(inv2);
    inv2_mat = reshape(inv2_mat,[3*num_dimensions,size(inv2,2)])';
    coord = inv2_mat(:,1:2*num_dimensions);
% for i=1:size(inv2,2)
%   coord = inv2{i}(1:2*num_dimensions);
%   filename = sprintf(out_format, i);
    filename = out_format;

  write_octave_like_output(filename,coord);
% end

out_format = strcat(element_dir, "element_table_3_4_FG.dat");
num_dimensions = 3;
lev_vec = zeros(num_dimensions,1) + 4;
grid_type = 'FG';
[fwd3, inv3] = create_hash_table(lev_vec, grid_type);
    inv3_mat = cell2mat(inv3);
    inv3_mat = reshape(inv3_mat,[3*num_dimensions,size(inv3,2)])';
    coord = inv3_mat(:,1:2*num_dimensions);
% for i=1:size(inv3,2)
%   coord = inv3{i}(1:2*num_dimensions);
%   filename = sprintf(out_format, i);
    filename = out_format;  
    write_octave_like_output(filename,coord);
% end
