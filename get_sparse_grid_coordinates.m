function [coordinates] = get_sparse_grid_coordinates(pde, opts, hash_table)
%%
% coordinates is a (num_pts x num_dims) array of real space coordinates of
% the sparse-grid element locations.

if opts.use_oldhash
    num_elems = numel(hash_table);
else
    num_elems = numel(hash_table.elements_idx);
end
num_dims  = numel(pde.dimensions);

%%
% Get center coordinates for each element

coordinates = zeros(num_elems,num_dims);
depth       = zeros (num_elems, 1);

plot_grid = false;
if plot_grid && numel(pde.dimensions)==2
    num_lev1 = max(hash_table.elements.lev_p1(:,1));
    num_lev2 = max(hash_table.elements.lev_p1(:,2));
end

z=0;
for elem=1:num_elems
    idx = hash_table.elements_idx(elem);
    coordinates(elem,:) = get_my_realspace_coord(pde, opts, hash_table, idx);
    
    if plot_grid && numel(pde.dimensions)==2
        depth(elem) = sum (hash_table.elements.lev_p1(idx,:)-1,'all');
        lev1 = hash_table.elements.lev_p1(idx,1)-1;
        lev2 = hash_table.elements.lev_p1(idx,2)-1;
        
        x_min = pde.dimensions{1}.min;
        x_max = pde.dimensions{1}.max;
        x_rng = x_max-x_min;
        
        y_min = pde.dimensions{2}.min;
        y_max = pde.dimensions{2}.max;
        y_rng = y_max-y_min;
        
        x = coordinates(elem,1);
        y = coordinates(elem,2);
        dx = 1/(2^max([lev1,1]))*x_rng;
        dy = 1/(2^max([lev2,1]))*y_rng;
        
        x1 = x+dx;
        x2 = x+dx;
        x3 = x-dx;
        x4 = x-dx;
        x5 = x1;
        
        y1 = y+dy;
        y2 = y-dy;
        y3 = y-dy;
        y4 = y+dy;
        y5 = y1;
        
        pos = num_lev1*lev2+lev1+1;
        subplot(num_lev1,num_lev2,pos)
        hold on
        plot([x1,x2,x3,x4,x5],[y1,y2,y3,y4,y5],'k','LineWidth',2,'MarkerSize',10);
        plot(x,y,'+k');
        xlim([x_min x_max])
        ylim([y_min y_max])
        lev_str = ['(',num2str(lev1),',',num2str(lev2),')'];
        set(gca,'FontSize',16)
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
        title(lev_str);
        
        assert(depth(elem)>=0);
        z = z - 1;
        
    end
    
end

if plot_grid
    figure
    plot(coordinates(:,1),coordinates(:,2),'+k','LineWidth',2)
    xlabel('p');
    ylabel('\zeta');
    lev_str = ['Lev = (',num2str(pde.dimensions{1}.lev),',',num2str(pde.dimensions{2}.lev),')'];
    title(lev_str);
    set(gca,'FontSize',16)
end

%%
% Get deg many points in each dimension for each element

% TODO


end
