function write_octave_like_output (file_name, mat)

%%
% Create a header string from present git HEAD hash

[status,git_hash] = system('git rev-parse HEAD');
s1 = '# Created from matlab git hash ';
s2 = git_hash;
header = [s1 s2];

%%
% Check to see if the file is changed by writing the file without the git
% hash header, running diff on it with the previous gold file, if there is
% a diff, then continue below.

root = get_root_folder();

gold_path = [root,'/gold/'];
gold_file_name = [gold_path,char(file_name)];
gold_file_name_tmp = [gold_file_name,'.tmp'];

update_file = false;

if exist(gold_file_name)
    
    existing_hash = get_existing_hash(gold_file_name);
    
    existing_header = [s1 existing_hash];
    
    write_file(gold_file_name_tmp, mat, existing_header);
    
    command = ['diff ', gold_file_name, ' ', gold_file_name_tmp];
    
    [status,cmdout] = system(command);
    
    delete(gold_file_name_tmp)
    
    if ~isempty(cmdout)
        update_file = true;
    end
    
else
    
    update_file = true;
    
end

% There is a diff or no file present
if update_file
    disp(['Updating gold data file ', gold_file_name]);
    write_file(gold_file_name, mat, header);
end

end

function write_file(file_name, mat, header)

write_header = false;
if nargin > 2
    write_header = true;
end

[num_rows,num_cols] = size(mat);

fid = fopen(file_name,'w');
if write_header
    fprintf(fid, header);
end
fprintf(fid,'# name: foo\n');
if isscalar(mat)
    fprintf(fid,'# type: scalar\n');
else
    fprintf(fid,'# type: matrix\n');
    fprintf(fid,'# rows: %i\n', num_rows);
    fprintf(fid,'# columns: %i\n', num_cols);
end

fclose(fid);

%%
% Write data

dlmwrite(file_name, mat, '-append', 'delimiter', ' ','precision',16);

end
