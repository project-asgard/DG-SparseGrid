function fileLength = get_file_length(fd)
% extracts file length in bytes from a file opened by fopen
% seek must be positioned at beginning on entry

seek = ftell(fd);
% move to end
fseek(fd, 0, 1);
% read end position
fileLength = ftell(fd);
% move to previous position
fseek(fd, seek, -1);
end
