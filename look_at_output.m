[root_directory] = get_root_folder();

files = dir([root_directory,'/output/asgard-out*.mat']);

num_files = numel(files);

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
    x = nodes_fixed_grid{1};
    y = f_fixed_grid(end-2,:);
    legend_strings{f} = [files(f).name,', DoF = ', num2str(dof)];
    if f>1; hold on; end 
    if contains(files(f).name,'FG')
        p=semilogy(x,y,'LineWidth',10);%,'Color','blue');
        p.Color(4)=0.3;
    else
        p=semilogy(x,y,'--','LineWidth',2);%,'Color','black');
        p.Color(4)=1.0;
    end
end
hold off
legend(legend_strings)

for f=1:num_files
    load(files(f).name);
    files(f).name
    if iscell(fval_t)
        [f_fixed_grid,nodes_fixed_grid] = get_realspace_fixed_grid_solution(pde,opts,fval_t{end},hash_table);
    else
        [f_fixed_grid,nodes_fixed_grid] = get_realspace_fixed_grid_solution(pde,opts,fval_t(:,end));
    end
    subplot(1,2,2)
    if f>1; hold on; end    
    x = nodes_fixed_grid{2};
    y = f_fixed_grid(:,floor(numel(x)/10));
    if contains(files(f).name,'FG')
        p=semilogy(x,y,'LineWidth',10);%,'Color','blue');
        p.Color(4)=0.3;
    else
        p=semilogy(x,y,'--','LineWidth',2);%,'Color','black');
        p.Color(4)=1.0;
    end   
end
hold off
legend(legend_strings)
disp('')
