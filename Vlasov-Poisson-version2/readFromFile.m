function [vector] = readFromFile(path)
%reads double precision, binary formatted vector from path

fd = fopen(path,'r');
%8 bytes in a double
n = getFileLength(fd) / 8;
vector = fread(fd,n,'double');
fclose(fd);

end

