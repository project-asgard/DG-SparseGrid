function connect_2D_create_test_data ()

Lev = [3,4,5]
Dim = 2;


for l=1:numel(Lev)
    
    [HASH,HASHInv] = HashTable(Lev(l),Dim);
    nHash = numel(HASHInv);
    
    Con2D = Connect2D(Lev(l),HASH,HASHInv);
    ConFlat = cell2mat(Con2D);
    
    levStr = sprintf("%1.1i",Lev(l));
    fname = 'tests/connect_2D/matlab-outputs-'+levStr+'.dat';
    
    fd = fopen(fname,'w'); % where file.dat is the name you want to save to
    fwrite(fd,ConFlat,'double'); % where U is the vector/matrix you want to store, double is the typename
    fclose(fd);
    
end


end
