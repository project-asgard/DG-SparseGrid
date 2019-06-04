function plot_adapt_triangle(pde,plot_num)

num_dimensions  = numel (pde.dimensions);
num_elements    = numel (pde.elementsIDX);

subplot(3,3,plot_num);
cla
hold on
if num_dimensions == 2
    
    for n=1:num_elements
        idx = pde.elementsIDX(n);
        lev_vec = pde.elements.lev_p1(idx,:)-1;
        plot(lev_vec(1),lev_vec(2),'o','MarkerSize',8,'MarkerFaceColor','k');
        
    end
end
hold off

end