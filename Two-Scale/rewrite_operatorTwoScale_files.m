function stat = rewrite_operatorTwoScale_files()

% List all mat files

s = what();

N = numel(s.mat);

options = '-ascii';

for f=1:N
    
    load(s.mat{f});
    
    str = extractAfter(s.mat{f},13);
    
    saveStrH0 = ['single/two_scale_rel_H0' str];
    saveStrG0 = ['single/two_scale_rel_G0' str];
    
    save(saveStrH0,'H0',options);
    save(saveStrG0,'G0',options);
    
end

end