function [] = write_A_Data_to_file(A_Data,Lev,Deg)

vMassV = A_Data.vMassV;
GradX  = A_Data.GradX;
GradV  = A_Data.GradV;
EMassX = A_Data.EMassX;

element_global_row_index = A_Data.element_global_row_index;
element_local_1_index = A_Data.element_local_1_index;
element_local_2_index = A_Data.element_local_2_index;
element_n_connected = A_Data.element_n_connected;

connected_global_col_index = A_Data.connected_global_col_index;
connected_local_1_index = A_Data.connected_local_1_index;
connected_local_2_index = A_Data.connected_local_2_index;

LevDegStr = sprintf('lev-%i_deg-%i',Lev,Deg);

if ~exist(['data'],'dir'),
   [status, msg, msgid] = mkdir('data');
   isok = (status == 1);
   if (~isok),
     error(sprintf('problem in creating directory %s,msg=%g,msgid=%g', ...
             'data', msg, msgid ));
   end;
end;

if ~exist(['data/' LevDegStr],'dir'), 
  [status, msg, msgid] = mkdir('data/',LevDegStr); 
  isok = (status == 1);
  if (~isok),
    error(sprintf('problem in creating director %s,msg=%g,msgid=%g', ...
            ['data/' LevDegStr], msg, msgid));
  end;
end;

options = '';
be_compatible_with_matlab = 0;
if (be_compatible_with_matlab),
  options = '-ascii';
end;

save(['data/' LevDegStr '/vMassV.mat'],'vMassV', options);
save(['data/' LevDegStr '/GradX.mat'],'GradX', options);
save(['data/' LevDegStr '/GradV.mat'],'GradV', options);
save(['data/' LevDegStr '/EMassX.mat'],'EMassX', options);

save(['data/' LevDegStr '/element_global_row_index.mat'],...
                         'element_global_row_index', options);
save(['data/' LevDegStr '/element_local_1_index.mat'],...
                        'element_local_1_index', options);
save(['data/' LevDegStr '/element_local_2_index.mat'],...
                        'element_local_2_index', options);
save(['data/' LevDegStr '/element_n_connected.mat'],...
                        'element_n_connected', options);

save(['data/' LevDegStr '/connected_global_col_index.mat'],...
                        'connected_global_col_index', options);
save(['data/' LevDegStr '/connected_local_1_index.mat'],...
                        'connected_local_1_index', options);
save(['data/' LevDegStr '/connected_local_2_index.mat'],...
                        'connected_local_2_index', options);

end
