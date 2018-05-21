function stat = OperatorTwoScale_create_test_data()

%% Run the OperatorTwoScale function to produce a set of input output pairs over which we can test the C++ implemetation. 

degList = [2,3,4,2];
levList = [3,4,5,6];

for i=1:numel(degList)
   
    output1 = OperatorTwoScale(degList(i),levList(i));
    
    fname = sprintf('test/fmwt_comp_%2.2i_%2.2i.mat',degList(i),levList(i));
    
    save(fname,'output1','-ascii'); 
    
end

end