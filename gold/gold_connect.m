%% Generate C++ testing data

% connectivity testing files
data_dir = strcat("generated-inputs", "/", "connectivity", "/");
root = get_root_folder();
[stat,msg] = mkdir ([root,'/gold/',char(data_dir)]);

% 1d indexing
out_format = strcat(data_dir, "get_1d_%d_%d.dat");
levs = [0, 0, 5];
cells = [0, 1, 9];
for i=1:size(levs,2)
index = lev_cell_to_1D_index(levs(i), cells(i));
filename = sprintf(out_format, levs(i), cells(i));
write_octave_like_output(filename,index);
end

% 1d connectivity
out_format = strcat(data_dir, "connect_1_%d.dat");
levs = [2, 3, 8];
for i=1:size(levs,2)
connectivity = full(connect_1D(levs(i)));
filename = sprintf(out_format, levs(i));
write_octave_like_output(filename,connectivity);
end

% nd connectivity
out_format = strcat(data_dir, "connect_n_2_3_FG_%d.dat");
dims = 2;
levs = 3;
grid = 'FG';
lev_sum = 6;
lev_max = 3;
lev_vec = zeros(dims,1)+levs;
[fwd, rev] = hash_table_nD(lev_vec, grid);
connectivity = connect_nD(dims, fwd, rev, lev_sum, lev_max);
for i=1:size(connectivity, 2)
filename = sprintf(out_format, i);
element = connectivity{i};
write_octave_like_output(filename,element);
end

out_format = strcat(data_dir, "connect_n_3_4_SG_%d.dat");
dims = 3;
levs = 4;
grid = 'SG';
lev_sum = 4;
lev_max = 4;
lev_vec = zeros(dims,1)+levs;
[fwd, rev] = hash_table_nD(lev_vec, grid);
connectivity = connect_nD(dims, fwd, rev, lev_sum, lev_max);
for i=1:size(connectivity, 2)
filename = sprintf(out_format, i);
element = connectivity{i};
write_octave_like_output(filename,element);
end

