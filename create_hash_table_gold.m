% Generate testing data for C++ 

addpath(genpath(pwd));

% element_table testing files
element_dir = strcat(pwd, "/", "generated-inputs", "/", "element_table", "/");
[status,msg] = mkdir (element_dir);

out_format = strcat(element_dir, "element_table_1_1_SG_%d.dat");
num_dimensions = 1;
lev_vec = zeros(num_dimensions,1) + 1;
grid_type = 'SG';
[fwd1, inv1] = create_hash_table(lev_vec, grid_type);
for i=1:size(inv1,2)
  coord = inv1{i}(1:2*num_dimensions);
  filename = sprintf(out_format, i);
  write_octave_like_output(filename,coord);
end

out_format = strcat(element_dir, "element_table_2_3_SG_%d.dat");
num_dimensions = 2;
lev_vec = zeros(num_dimensions,1) + 3;
grid_type = 'SG';
[fwd2, inv2] = create_hash_table(lev_vec, grid_type);
for i=1:size(inv2,2)
  coord = inv2{i}(1:2*num_dimensions);
  filename = sprintf(out_format, i);
  write_octave_like_output(filename,coord);
end

out_format = strcat(element_dir, "element_table_3_4_FG_%d.dat");
num_dimensions = 3;
lev_vec = zeros(num_dimensions,1) + 4;
grid_type = 'FG';
[fwd3, inv3] = create_hash_table(lev_vec, grid_type);
for i=1:size(inv3,2)
  coord = inv3{i}(1:2*num_dimensions);
  filename = sprintf(out_format, i);
  write_octave_like_output(filename,coord);
end
