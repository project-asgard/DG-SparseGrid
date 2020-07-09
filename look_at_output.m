[root_directory] = get_root_folder();

files = dir([root_directory,'/output/SG-comparison-d4/asgard-out*.mat']);

num_files = numel(files);

figure('position',[0 0 1400 500])
subplot(1,2,1)

colors = ["black","blue","red","green","cyan","magenta","black","blue","red","green","cyan","magenta"];
fg_cnt = 1;
sg_cnt = 1;

% load correct solution
correct = [root_directory,'/output/','FG-convergence-l4/asgard-out-fokkerplanck2_complete-l4-d15-FG-10-5-BE-TIA-true-TIBA-false.mat'];
load(correct);
if iscell(fval_t)
    [f_fixed_grid,nodes_fixed_grid] = get_realspace_fixed_grid_solution(pde,opts,fval_t{end},hash_table);
    dof = numel(fval_t{end});
else
    [f_fixed_grid,nodes_fixed_grid] = get_realspace_fixed_grid_solution(pde,opts,fval_t(:,end));
    dof = numel(fval_t(:,end));
end
f_correct = f_fixed_grid;

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
        p=semilogy(x,y,'LineWidth',5,'Color',colors(fg_cnt));
        fg_cnt = fg_cnt + 1;
        p.Color(4)=0.3;
        error_fg(fg_cnt) = sqrt(mean((f_correct - f_fixed_grid).^2,'all'));
%         error_fg(fg_cnt) = norm(f_correct - f_fixed_grid);

        dof_fg(fg_cnt) = dof;
%         p.Color(4)=max([min([dof./50000 1]) 0.1]);
    else
        grid_str = 'SG';
        p=semilogy(x,y,'--','LineWidth',2,'Color',colors(sg_cnt));
        sg_cnt = sg_cnt+1;
        p.Color(4)=1.0;
        error_sg(sg_cnt) = sqrt(mean((f_correct - f_fixed_grid).^2,'all'));
%         error_sg(sg_cnt) = norm(f_correct - f_fixed_grid);

        dof_sg(sg_cnt) = dof;
    end
    xlabel('p');
    ylabel(['f(\zeta=',num2str(zeta),',p)']);
    legend_strings{f} = [grid_str,lev_str,', DoF = ', num2str(dof)];

end
hold off
legend(legend_strings,'FontSize',12)
set(gca,'FontSize',16)

fg_cnt = 1;
sg_cnt = 1;
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
        p=semilogy(x,y,'LineWidth',5,'Color',colors(fg_cnt));
        fg_cnt = fg_cnt + 1;
        p.Color(4)=0.3;
    else
        grid_str = 'SG';
        p=semilogy(x,y,'--','LineWidth',2,'Color',colors(sg_cnt));
        sg_cnt = sg_cnt+1;
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

% plot & save error

[dof_fg,ii] = sort(dof_fg);
error_fg = error_fg(ii);
[dof_sg,ii] = sort(dof_sg);
error_sg = error_sg(ii);

figure
loglog(dof_sg,error_sg,'-o','LineWidth',2)
xlabel('Degrees of Freedom');
ylabel('RMS Absolute Error');
hold on
loglog(dof_fg,error_fg,'-o','LineWidth',2)
hold off
set(gca,'FontSize',16)
dof_legend_strings{1} = ['Sparse Grid (deg=',num2str(pde.deg),')'];
dof_legend_strings{2} = ['Full Grid (deg=',num2str(pde.deg),')'];
legend(dof_legend_strings,'FontSize',12)


