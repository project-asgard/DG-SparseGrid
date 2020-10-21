function plot_adapt_triangle(pde,opts,hash_table,plot_num)

num_dimensions  = numel (pde.dimensions);
num_elements    = numel (hash_table.elements_idx);

if num_dimensions == 2
    subplot(4,3,plot_num);
    cla
    hold on
    
    for n=1:num_elements
        idx = hash_table.elements_idx(n);
        lev_vec = hash_table.elements.lev_p1(idx,:)-1;
        plot(lev_vec(1),lev_vec(2),'o','MarkerSize',8,'MarkerFaceColor','k');
        
    end
    hold off
    
end

end
