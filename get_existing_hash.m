function [hash] = get_existing_hash(file_name)

fid = fopen(file_name,'r');

[A,count] = fscanf(fid,'%s',7);

hash = [extractAfter(A,25) '\n'];

end
