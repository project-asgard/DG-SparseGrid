function write_octave_like_output (file_name, mat)

[num_rows,num_cols] = size(mat);

%%
% Write octave like header

[status,git_hash] = system('git rev-parse HEAD');
s1 = '# Created from matlab git hash ';
s2 = git_hash;
created_from_string = [s1 s2];

fid = fopen(file_name,'w');
fprintf(fid,created_from_string);
fprintf(fid,'# name: foo\n');
fprintf(fid,'# type: matrix\n');
fprintf(fid,'# rows: %i\n', num_rows);
fprintf(fid,'# columns: %i\n', num_cols);
fclose(fid);

%%
% Write data

dlmwrite(file_name, mat, '-append', 'delimiter', ' '); 

end