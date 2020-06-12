[root_directory] = get_root_folder();

files = dir([root_directory,'/output/asgard-out*.mat']);

num_files = numel(files);

subplot(1,2,1)

for f=1:num_files    
    load(files(f).name);
    x = nodes_fixed_grid_nodups{1};
    y = f_realspace_nD_fixed_grid_nodups_t(end-1,:,end);
    dof = numel(fval_t(:,1));
    legend_strings{f} = [files(f).name,', DoF = ', num2str(dof)];
    if f>1; hold on; end    
    p=semilogy(x,y,'LineWidth',5);%,'Color','black');
    p.Color(4)=0.5;
end
hold off
legend(legend_strings)

for f=1:num_files
    load(files(f).name);
    subplot(1,2,2)
    if f>1; hold on; end    
    x = nodes_fixed_grid_nodups{2};
    y = f_realspace_nD_fixed_grid_nodups_t(:,floor(numel(x)/10),end);
    p2=semilogy(x,y,'LineWidth',5);%,'Color','black');
    p2.Color(4)=0.5;    
end
hold off
legend(legend_strings)
disp('')
