function connect_1D_create_test_data ()

% Run in Vlasov-Poisson-version2 folder

Lev = [3,4,5];

for i=1:numel(Lev)

Con = full( Connect1D(Lev(i)) );

levStr = sprintf("%1.1i",Lev(i));
fname = 'tests/connect_1D/matlab-outputs-'+levStr+'.dat';

fd = fopen(fname,'w'); % where file.dat is the name you want to save to
fwrite(fd,Con,'double'); % where U is the vector/matrix you want to store, double is the typename
fclose(fd);

end


end