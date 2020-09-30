[root_directory] = get_root_folder();

files = dir([root_directory,'/output/asgard-out*.mat']);

num_files = numel(files);

for f=1:num_files
    
    load(files(f).name);
    x = nodes{1};
    y = f_realspace_nD_t(end,:,end);
    dof = numel(fval_t(:,1));
    legend_strings{f} = [files(f).name,', DoF = ', num2str(dof)];
    if f>1; hold on; end
    p=semilogy(x,y,'LineWidth',5);%,'Color','black');
    p.Color(4)=0.5;
    
end
hold off
legend(legend_strings)

