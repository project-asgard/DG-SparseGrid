%% Generate C++ testing data for matlab_utilities component

data_dir = strcat("generated-inputs", "/", "matlab_utilities", "/");
root = get_root_folder();
[stat,msg] = mkdir ([root,'/gold/',char(data_dir)]);

% these are used to test linspace()

out_format = strcat(data_dir, "linspace_");
w1 = linspace(-1,1,9);
w2 = linspace(1,-1,9);
w3 = linspace(-1,1,8);
write_octave_like_output(strcat(out_format,'neg1_1_9.dat'), w1);
write_octave_like_output(strcat(out_format,'1_neg1_9.dat'), w2);
write_octave_like_output(strcat(out_format,'neg1_1_8.dat'), w3);

% these are used to test read_vector_from_bin_file()

out_format = strcat(data_dir, "read_vector_bin_");
w = linspace(-1,1);
wT = w';
write_to_file(strcat(out_format, 'neg1_1_100.dat'), w);
write_to_file(strcat(out_format, 'neg1_1_100T.dat'), wT);

% these are used to test read_vector_from_txt_file()

out_format = strcat(data_dir, "read_vector_txt_");
w2 = linspace(-1,1);
w2T = w';
write_octave_like_output(strcat(out_format, 'neg1_1_100.dat'), w2);
write_octave_like_output(strcat(out_format, 'neg1_1_100T.dat'), w2T);

% these are used to test read_vector_from_txt_file()

for i = 0:4; for j = 0:4
  m(i+1,j+1) = 17/(i+1+j);
end
end

write_octave_like_output(strcat(data_dir, 'read_matrix_txt_5x5.dat'), m);

% test read_scalar_from_txt_file()

s = 42;
write_octave_like_output(strcat(data_dir, 'read_scalar_42.dat'), s);
