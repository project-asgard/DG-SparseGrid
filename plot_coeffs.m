function [res] = plot_coeffs(num_dims,max_lev,hash_table,deg,fval,num_rows,num_cols,pos,...
    refine_threshold, coarsen_threshold)

element_DOF = deg^num_dims;
num_elements = numel(fval)/element_DOF;
hold off
subplot(num_rows,num_cols,pos)
fval_element = zeros(num_elements,1);
depth = zeros(num_elements,1);
fval_depth = zeros(numel(fval),1);
for i=1:num_elements
    view = fval((i-1).*element_DOF+1:i*element_DOF);
    fval_element(i) = sqrt(sum(view.^2));
    lev_vec = md_idx_to_lev_pos(num_dims,max_lev,hash_table.elements_idx(i));
    depth(i) = sum(lev_vec);
    fval_depth((i-1).*element_DOF+1:i*element_DOF) = depth(i);
end

semilogy(fval_depth,abs(fval),'o','Color','black');
hold on
semilogy(depth,fval_element,'o','Color','blue');
yline(refine_threshold);
yline(coarsen_threshold);
hold off

end
