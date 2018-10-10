function stat = rewrite_operatorTwoScale_files()

% Run in the "Two-Scale" directory

% Generate all the original two_scale_rel_deg.mat files up to this deg ...

deg = 100;

options = '-ascii';

for k=1:deg
    
    [H0,G0,scale_co,phi_co]=MultiwaveletGen(k);
    
    saveStr = ['two_scale_rel_' num2str(k) '.mat'];
    save(saveStr,'H0','G0','scale_co','phi_co',options);
    
    fNameH0 = ['single/two_scale_rel_H0_' num2str(k) '.dat'];
    fNameG0 = ['single/two_scale_rel_G0_' num2str(k) '.dat'];
    fName_phi_co = ['single/two_scale_rel_phi_co_' num2str(k) '.dat'];
    fName_scale_co = ['single/two_scale_rel_scale_co_' num2str(k) '.dat'];
    
    fd = fopen(fNameH0,'w'); % where file.dat is the name you want to save to
    fwrite(fd,H0,'double'); % where U is the vector/matrix you want to store, double is the typename
    fclose(fd);
    
    fd = fopen(fNameG0,'w'); % where file.dat is the name you want to save to
    fwrite(fd,G0,'double'); % where U is the vector/matrix you want to store, double is the typename
    fclose(fd);    
    
    fd = fopen(fName_phi_co,'w'); % where file.dat is the name you want to save to
    fwrite(fd,phi_co,'double'); % where U is the vector/matrix you want to store, double is the typename
    fclose(fd);   
    
    fd = fopen(fName_scale_co,'w'); % where file.dat is the name you want to save to
    fwrite(fd,scale_co,'double'); % where U is the vector/matrix you want to store, double is the typename
    fclose(fd); 
end


end