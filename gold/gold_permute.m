% Generate testing data for C++ 

% permutations testing files

data_dir = strcat("generated-inputs", "/", "permutations", "/");
root = get_root_folder();
[stat,msg] = mkdir ([root,'/gold/',char(data_dir)]);

dims = [1, 2, 4, 6];
ns = [1, 4, 6, 8];
ords = [0, 1, 0, 1];

% perm leq
out_format = strcat(data_dir, "perm_leq_%d_%d_%d.dat");
count_out_format = strcat(data_dir, "perm_leq_%d_%d_%d_count.dat");
for i=1:size(dims,2)
  tuples = perm_leq(dims(i), ns(i), ords(i));
  count = perm_leq_count(dims(i), ns(i));
  filename = sprintf(out_format, dims(i), ns(i), ords(i));
  count_filename = sprintf(count_out_format, dims(i), ns(i), ords(i));
  write_octave_like_output(filename,tuples);
  write_octave_like_output(count_filename,count);
end

%perm eq
out_format = strcat(data_dir, "perm_eq_%d_%d_%d.dat");
count_out_format = strcat(data_dir, "perm_eq_%d_%d_%d_count.dat");
for i=1:size(dims,2)
  tuples = perm_eq(dims(i), ns(i), ords(i));
  count = perm_eq_count(dims(i), ns(i));
  filename = sprintf(out_format, dims(i), ns(i), ords(i));
  count_filename = sprintf(count_out_format, dims(i), ns(i), ords(i));
  write_octave_like_output(filename,tuples);
  write_octave_like_output(count_filename,count);
end

%perm max
out_format = strcat(data_dir, "perm_max_%d_%d_%d.dat");
count_out_format = strcat(data_dir, "perm_max_%d_%d_%d_count.dat");
for i=1:size(dims,2)
  tuples = perm_max(dims(i), ns(i), ords(i));
  count = perm_max_count(dims(i), ns(i));
  filename = sprintf(out_format, dims(i), ns(i), ords(i));
  count_filename = sprintf(count_out_format, dims(i), ns(i), ords(i));
  write_octave_like_output(filename,tuples);
  write_octave_like_output(count_filename,count);
end

%index_leq_max

filename = strcat(data_dir, "index_leq_max_4d_10s_4m.dat");
count_filename = strcat(data_dir, "index_leq_max_4d_10s_4m_count.dat");

level_sum = 10;
level_max = 4;
dim = 4;
lists{1} = 2:3;
lists{2} = 0:4;
lists{3} = 0:3;
lists{4} = 1:5;
result = index_leq_max(dim, lists, level_sum, level_max);
count = index_leq_max_count(dim, lists, level_sum, level_max);
write_octave_like_output(filename, result);
write_octave_like_output(count_filename, count);