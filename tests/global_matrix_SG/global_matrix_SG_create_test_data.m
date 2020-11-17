function stat = global_matrix_SG_create_test_data ()

Lev = 3;
Dim = 2;
Deg = 2;
compression = 4;

[HASH,HASHInv] = HashTable(Lev,Dim);

Con2D = Connect2D(Lev,HASH,HASHInv);

A_data = GlobalMatrixSG_SlowVersion(HASHInv,Con2D,Deg,compression);

fname = 'tests/global_matrix_SG/matlab-outputs-adata1.dat';
fd = fopen(fname,'w'); % where file.dat is the name you want to save to
fwrite(fd,A_data.element_global_row_index,'int'); % where U is the vector/matrix you want to store, double is the typename
fclose(fd);

fname = 'tests/global_matrix_SG/matlab-outputs-adata2.dat';
fd = fopen(fname,'w'); % where file.dat is the name you want to save to
fwrite(fd,A_data.element_local_1_index,'int'); % where U is the vector/matrix you want to store, double is the typename
fclose(fd);

fname = 'tests/global_matrix_SG/matlab-outputs-adata3.dat';
fd = fopen(fname,'w'); % where file.dat is the name you want to save to
fwrite(fd,A_data.element_local_2_index,'int'); % where U is the vector/matrix you want to store, double is the typename
fclose(fd);

fname = 'tests/global_matrix_SG/matlab-outputs-adata4.dat';
fd = fopen(fname,'w'); % where file.dat is the name you want to save to
fwrite(fd,A_data.element_n_connected,'int'); % where U is the vector/matrix you want to store, double is the typename
fclose(fd);

fname = 'tests/global_matrix_SG/matlab-outputs-adata5.dat';
fd = fopen(fname,'w'); % where file.dat is the name you want to save to
fwrite(fd,A_data.connected_global_col_index,'int'); % where U is the vector/matrix you want to store, double is the typename
fclose(fd);

fname = 'tests/global_matrix_SG/matlab-outputs-adata6.dat';
fd = fopen(fname,'w'); % where file.dat is the name you want to save to
fwrite(fd,A_data.connected_local_1_index,'int'); % where U is the vector/matrix you want to store, double is the typename
fclose(fd);

fname = 'tests/global_matrix_SG/matlab-outputs-adata7.dat';
fd = fopen(fname,'w'); % where file.dat is the name you want to save to
fwrite(fd,A_data.connected_local_2_index,'int'); % where U is the vector/matrix you want to store, double is the typename
fclose(fd);

end
