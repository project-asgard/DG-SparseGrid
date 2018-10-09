function stat = rewrite_operatorTwoScale_files()

% % Generate all the original two_scale_rel_deg.mat files up to this deg ...
% 
% deg = 100;
% 
% for k=1:deg
%     
%     [H0,G0,scale_co,phi_co]=MultiwaveletGen(k);
%     
%     saveStr = ['two_scale_rel_',num2str(k)];
%     save(saveStr,'H0','G0','scale_co','phi_co');
%     
% end

% List all mat files

s = what();

N = numel(s.mat);

options = '-ascii';

for f=1:N
    
    load(s.mat{f});
    
    str = extractAfter(s.mat{f},13);
    
    saveStrH0 = ['single/two_scale_rel_H0' str];
    saveStrG0 = ['single/two_scale_rel_G0' str];
    saveStr_phi_co = ['single/two_scale_rel_phi_co' str];
    
    save(saveStrH0,'H0',options);
    save(saveStrG0,'G0',options);
    save(saveStr_phi_co,'phi_co',options);
    
end

end