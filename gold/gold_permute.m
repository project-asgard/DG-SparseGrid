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


% perm_leq_d / perm_eq_d

out_format = strcat(data_dir, "perm_leq_d_%d_%d.dat");
count_out_format = strcat(data_dir, "perm_leq_d_%d_%d_count.dat");
e_out_format = strcat(data_dir, "perm_eq_d_%d_%d.dat");
e_count_out_format = strcat(data_dir, "perm_eq_d_%d_%d_count.dat");

levels{1} = [3, 3];
levels{2} = [1, 4];
levels{3} = [1, 5, 8];
levels{4} = [10, 6, 9, 10];
levels{5} = [2, 10, 1, 5, 4, 7];

for i=1:size(levels, 2)
   sort = mod(i, 2);

   count = perm_leq_d_count(size(levels{i}, 2),levels{i}, max(levels{i}));
   result = perm_leq_d(size(levels{i}, 2),levels{i},max(levels{i}),sort);
   
   disp(levels{i});
   disp(size(levels{i}, 2));
   disp(sum(levels{i}));
   e_count = perm_eq_d_count(size(levels{i}, 2), levels{i}, max(levels{i}));
   e_result = perm_eq_d(size(levels{i}, 2),levels{i}, max(levels{i}),sort);
   
   filename = sprintf(out_format, size(levels{i}, 2), sort);
   count_filename = sprintf(count_out_format, size(levels{i}, 2), sort);   
   write_octave_like_output(filename,result);
   write_octave_like_output(count_filename,count);
   
   filename = sprintf(e_out_format, size(levels{i}, 2), sort);
   count_filename = sprintf(e_count_out_format, size(levels{i}, 2), sort);   
   write_octave_like_output(filename,e_result);
   write_octave_like_output(count_filename,e_count);
end




