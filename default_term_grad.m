term_grad.type = 'grad';
term_grad.G = @(x,p,t,dat) x*0+1;
term_grad.TD = 0;
term_grad.dat = [];
term_grad.LF = 0; % central flux
term_grad.name = 'grad';
term_grad.BCL = 'N';
term_grad.BCR = 'N';
for d=1:num_dimensions % BC variation in all dimensions
    term_grad.BCL_fList{d} = @(x,p,t) x.*0;
    term_grad.BCR_fList{d} = @(x,p,t) x.*0;
end
term_grad.BCL_fList{num_dimensions+1} = @(t,p) 1;  % time variation
term_grad.BCR_fList{num_dimensions+1} = @(t,p) 1;
