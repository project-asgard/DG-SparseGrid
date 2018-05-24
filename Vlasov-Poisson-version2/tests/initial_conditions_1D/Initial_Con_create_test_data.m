function stat = Initial_Con_create_test_data ()

%% Write input and output test data for Initial_Con.m using the Vlasov4 case.
% 
% Run from Vlasov-Poission-version2 folder

Lev_x = 3;
Lev_v = 3;
Lmax = 20.944;
Vmax = 13.0;
k = 2;

fname = 'tests/initial_conditions_1D/matlab-inputs.mat';
load(fname);

[f_x,f_v] = Intial_Con(Lev_x,Lev_v,k,Lmax,Vmax,PDE,FMWT_COMP_x,FMWT_COMP_v);

fname = 'tests/initial_conditions_1D/matlab-inputs-FMWTx.dat';

fd = fopen(fname,'w'); % where file.dat is the name you want to save to
fwrite(fd,FMWT_COMP_x,'double'); % where U is the vector/matrix you want to store, double is the typename
fclose(fd);

fname = 'tests/initial_conditions_1D/matlab-inputs-FMWTv.dat';

fd = fopen(fname,'w'); % where file.dat is the name you want to save to
fwrite(fd,FMWT_COMP_v,'double'); % where U is the vector/matrix you want to store, double is the typename
fclose(fd);


fname = 'tests/initial_conditions_1D/matlab-outputs-f_x.dat';

fd = fopen(fname,'w'); % where file.dat is the name you want to save to
fwrite(fd,f_x,'double'); % where U is the vector/matrix you want to store, double is the typename
fclose(fd);

fname = 'tests/initial_conditions_1D/matlab-outputs-f_v.dat';

fd = fopen(fname,'w'); % where file.dat is the name you want to save to
fwrite(fd,f_v,'double'); % where U is the vector/matrix you want to store, double is the typename
fclose(fd);

end