function [] = write_A_data_to_file(A_data,Lev,Deg)

vMassV = A_data{1}.vMassV;
GradX  = A_data{1}.GradX;
GradV  = A_data{1}.GradV;
EMassX = A_data{1}.EMassX;

element_global_row_index = A_data{1}.element_global_row_index;
element_local_1_index = A_data{1}.element_local_1_index;
element_local_2_index = A_data{1}.element_local_2_index;
element_n_connected = A_data{1}.element_n_connected;

connected_global_col_index = A_data{1}.connected_global_col_index;
connected_local_1_index = A_data{1}.connected_local_1_index;
connected_local_2_index = A_data{1}.connected_local_2_index;

LevDegStr = sprintf('lev-%i_deg-%i',Lev,Deg);

if ~exist(['data/' LevDegStr]); mkdir('data/',LevDegStr); end

save(['data/' LevDegStr '/vMassV.mat'],'vMassV');
save(['data/' LevDegStr '/GradX.mat'],'GradX');
save(['data/' LevDegStr '/GradV.mat'],'GradV');
save(['data/' LevDegStr '/EMassX.mat'],'EMassX');

save(['data/' LevDegStr '/element_global_row_index.mat'],'element_global_row_index');
save(['data/' LevDegStr '/element_local_1_index.mat'],'element_local_1_index');
save(['data/' LevDegStr '/element_local_2_index.mat'],'element_local_2_index');
save(['data/' LevDegStr '/element_n_connected.mat'],'element_n_connected');

save(['data/' LevDegStr '/connected_global_col_index.mat'],'connected_global_col_index');
save(['data/' LevDegStr '/connected_local_1_index.mat'],'connected_local_1_index');
save(['data/' LevDegStr '/connected_local_2_index.mat'],'connected_local_2_index');

end