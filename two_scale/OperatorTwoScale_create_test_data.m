function stat = OperatorTwoScale_create_test_data()

%% Run the OperatorTwoScale function to produce a set of input output pairs over which we can test the C++ implemetation.

degList = [2,3,4,2];
levList = [3,4,5,6];

for i=1:numel(degList)
    
    output = OperatorTwoScale(degList(i),levList(i));
    
    fname = sprintf('test/fmwt_comp_%2.2i_%2.2i.dat',degList(i),levList(i));
    
    %save(fname,'output','-ascii');
    
    % Write files compatible with the ReadVectorFromFile function in
    % matrixIO.hpp to read this output (file.dat) into a C++ vector
    
    fd = fopen(fname,'w'); % where file.dat is the name you want to save to
    fwrite(fd,output,'double'); % where U is the vector/matrix you want to store, double is the typename
    fclose(fd);
    
end

end
