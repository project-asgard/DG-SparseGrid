[root_directory] = get_root_folder();

files = dir([root_directory,'/output/asgard-out*.mat']);

num_files = numel(files);

figure('position',[0 0 1400 500])
subplot(1,2,1)

for f=1:num_files    
    load(files(f).name);
    files(f).name
    if iscell(fval_t)
        [f_fixed_grid,nodes_fixed_grid] = get_realspace_fixed_grid_solution(pde,opts,fval_t{end},hash_table);
        dof = numel(fval_t{end});
    else
        [f_fixed_grid,nodes_fixed_grid] = get_realspace_fixed_grid_solution(pde,opts,fval_t(:,end));
        dof = numel(fval_t(:,end));
    end
    zeta = nodes_fixed_grid{2}(end-2);
    x = nodes_fixed_grid{1};
    y = f_fixed_grid(end-2,:);
    lev_str = [', lev=',num2str(pde.dimensions{1}.lev),' deg=',num2str(pde.deg)];
    if f>1; hold on; end 
    if contains(files(f).name,'FG')
        grid_str = 'FG';
        p=semilogy(x,y,'LineWidth',5,'Color','black');
        p.Color(4)=min([dof./50000 1]);
    else
        grid_str = 'SG';
        p=semilogy(x,y,'--','LineWidth',2);%,'Color','black');
        p.Color(4)=1.0;
    end
    xlabel('p');
    ylabel(['f(\zeta=',num2str(zeta),',p)']);
    legend_strings{f} = [grid_str,lev_str,', DoF = ', num2str(dof)];

end
hold off
legend(legend_strings,'FontSize',12)
set(gca,'FontSize',16)

for f=1:num_files
    load(files(f).name);
    files(f).name
    if iscell(fval_t)
        [f_fixed_grid,nodes_fixed_grid] = get_realspace_fixed_grid_solution(pde,opts,fval_t{end},hash_table);
        dof = numel(fval_t{end});
    else
        [f_fixed_grid,nodes_fixed_grid] = get_realspace_fixed_grid_solution(pde,opts,fval_t(:,end));
        dof = numel(fval_t(:,end));
    end
    subplot(1,2,2)
    if f>1; hold on; end  
    pcoord = floor(numel(x)/4);
    pval = nodes_fixed_grid{1}(pcoord);
    x = nodes_fixed_grid{2};
    y = f_fixed_grid(:,floor(pcoord));
    y(y<=0)=1e-12;
    if contains(files(f).name,'FG')
        grid_str = 'FG';
        p=semilogy(x,y,'LineWidth',5,'Color','black');
        p.Color(4)=min([dof./50000 1]);
    else
        grid_str = 'SG';
        p=semilogy(x,y,'--','LineWidth',2);%,'Color','black');
        p.Color(4)=1.0;
    end  
    xlabel('\zeta');
    ylabel(['f(\zeta',',p=',num2str(pval),')']);
end
hold off
legend(legend_strings,'FontSize',12)
ylim([1e-10 1])
set(gca,'FontSize',16)
exportgraphics(gcf,'asgard-2d-slice.pdf','ContentType','vector','BackgroundColor','none');
disp('')
