[root_directory] = get_root_folder();

% files = dir([root_directory,'/output/SG-comparison-d4/asgard-out*.mat']);
% files = dir([root_directory,'/output/FG-convergence-d3/asgard-out*.mat']);
% files = dir([root_directory,'/output/asgard-out*.mat']);
% files = dir([root_directory,'/output/SG-comparison-d4/asgard-out*.mat']);
% files = dir([root_directory,'/output/FG-convergence-l4/asgard-out*.mat']);
files = dir([root_directory,'/output/ASG-comparison-d4-max/asgard-out*.mat']);
% files = dir([root_directory,'/output/ASG-comparison-d6/asgard-out*.mat']);

num_files = numel(files);

figure('position',[0 0 1400 500])

colors = ["black","blue","red","green","cyan","magenta","black","blue","red","green","cyan","magenta"];
fg_cnt = 1;
sg_cnt = 1;
fg_f = {};
sg_f = {};
dof_fg = [];
dof_sg = [];

p_err_idx = 200;
z_err_idx = 254;

% load correct solution
correct = [root_directory,'/output/','FG-convergence-l4/asgard-out-fokkerplanck2_complete-l4-d12-FG-10-5-BE-TIA-true-TIBA-false.mat'];
% correct = [root_directory,'/output/','asgard-out-l4-d4-SG-dt10-adapt-1.0e-06.mat'];
% correct = [root_directory,'/output/','ASG-comparison-d6/stash/asgard-out-l6-d6-SG-dt10-adapt-1.0e-07-max.mat'];
clear pde; clear opts;
load(correct);
if isfield(pde,'deg') % catch for older output files
    opts.deg = pde.deg;
    opts.lev_vec = pde.lev_vec;
    opts.max_lev = pde.max_lev;
end

if iscell(fval_t)
    [f_fixed_grid,nodes_fixed_grid] = get_realspace_fixed_grid_solution(pde,opts,fval_t{end},hash_table);
    dof = numel(fval_t{end});
else
    [f_fixed_grid,nodes_fixed_grid] = get_realspace_fixed_grid_solution(pde,opts,fval_t(:,end));
    dof = numel(fval_t(:,end));
end
f_correct = f_fixed_grid;

z_idx = numel(nodes_fixed_grid{2})-2;
p_idx = floor(numel(nodes_fixed_grid{1})/4);

subplot(1,2,1)
x = nodes_fixed_grid{1};
y = f_correct(z_idx,:);
p=semilogy(x,y,'LineWidth',10);
p.Color(4)=0.5;

subplot(1,2,2)
x = nodes_fixed_grid{2};
y = f_correct(:,floor(p_idx));
y(y<=0)=1e-12;
p=semilogy(x,y,'LineWidth',10);
p.Color(4)=0.5;

legend_strings{1} = ['correct'];


% load other solutions

for f=1:num_files  
    clear pde; clear opts;
    load(files(f).name);
    if isfield(pde,'deg') % catch for older output files
        opts.deg = pde.deg;
        opts.lev_vec = pde.lev_vec;
        opts.max_lev = pde.max_lev;
    end
    files(f).name
    if iscell(fval_t)
        [f_fixed_grid,nodes_fixed_grid] = get_realspace_fixed_grid_solution(pde,opts,fval_t{end},hash_table);
        dof = numel(fval_t{end});
    else
        [f_fixed_grid,nodes_fixed_grid] = get_realspace_fixed_grid_solution(pde,opts,fval_t(:,end));
        dof = numel(fval_t(:,end));
    end

    zeta = nodes_fixed_grid{2}(z_idx);
    x = nodes_fixed_grid{1};
    y = f_fixed_grid(z_idx,:);
    lev_str = [', lev=',num2str(pde.dimensions{1}.lev),' deg=',num2str(opts.deg)];
    subplot(1,2,1)
    hold on 
    if contains(files(f).name,'FG')
        grid_str = 'FG';
        p=semilogy(x,y,'LineWidth',5,'Color',colors(fg_cnt));
%         error_fg(fg_cnt) = sqrt(mean((f_correct - f_fixed_grid).^2,'all'));
        error_fg(fg_cnt) = norm(f_correct - f_fixed_grid,Inf);
        dof_fg(fg_cnt) = dof; 
        fg_f{fg_cnt} = f_fixed_grid;
        fg_files{fg_cnt} = files(f).name;
        deg_fg(fg_cnt) = opts.deg;

        fg_cnt = fg_cnt + 1;
        p.Color(4)=0.3;
%         p.Color(4)=max([min([dof./50000 1]) 0.1]);
    else
        grid_str = 'SG';
        p=semilogy(x,y,'--','LineWidth',2,'Color',colors(sg_cnt));
        p.Color(4)=1.0;
%         error_sg(sg_cnt) = sqrt(mean((f_correct - f_fixed_grid).^2,'all'));
%         error_sg(sg_cnt) = sqrt(mean((f_correct(z_err_idx,p_err_idx) - f_fixed_grid(z_err_idx,p_err_idx)).^2,'all'));
        error_sg(sg_cnt) = norm(f_correct - f_fixed_grid,Inf);
        
%         levels = [10,1,0.1,0.01,0.001,0.0001,0.00001,0.000001,0.0000001];
%         figure(33)
%         subplot(1,3,1)
%         contourf(real(log10(f_correct)));
%         colorbar()
%         set(gca,'fontsize', 18);
%         
%         subplot(1,3,2)
%         contourf(real(log10(f_fixed_grid)));
%         colorbar()
%         set(gca,'fontsize', 18);
%         
%         subplot(1,3,3)
%         contourf(f_correct-f_fixed_grid);
%         colorbar()
%         set(gca,'fontsize', 18);
      
        dof_sg(sg_cnt) = dof;
        sg_f{sg_cnt} = f_fixed_grid;
        sg_files{sg_cnt} = files(f).name;
        sg_cnt = sg_cnt+1;
    end
    xlabel('p');
    ylabel(['f(\zeta=',num2str(zeta),',p)']);
    legend_strings{f+1} = [grid_str,lev_str,', DoF = ', num2str(dof)];

end
hold off
legend(legend_strings,'FontSize',12)
set(gca,'FontSize',16)

fg_cnt = 1;
sg_cnt = 1;
for f=1:num_files
    clear pde; clear opts;
    load(files(f).name);
    if isfield(pde,'deg') % catch for older output files
        opts.deg = pde.deg;
        opts.lev_vec = pde.lev_vec;
        opts.max_lev = pde.max_lev;
    end
    files(f).name
    if iscell(fval_t)
        [f_fixed_grid,nodes_fixed_grid] = get_realspace_fixed_grid_solution(pde,opts,fval_t{end},hash_table);
        dof = numel(fval_t{end});
    else
        [f_fixed_grid,nodes_fixed_grid] = get_realspace_fixed_grid_solution(pde,opts,fval_t(:,end));
        dof = numel(fval_t(:,end));
    end
    subplot(1,2,2)
    hold on  
    pval = nodes_fixed_grid{1}(p_idx);
    x = nodes_fixed_grid{2};
    y = f_fixed_grid(:,floor(p_idx));
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

if ~isempty(dof_fg)
    [dof_fg,ii] = sort(dof_fg);
    error_fg = error_fg(ii);
    fg_files = fg_files(ii);
    fg_f = fg_f(ii);
    deg_fg = deg_fg(ii);
    fg_convergence = [];
    if numel(fg_f)>1
        for i=2:numel(fg_f)
            fg_diff = fg_f{i}-fg_f{i-1};
%             fg_diff = fg_f{i}(z_err_idx,p_err_idx)-fg_f{i-1}(z_err_idx,p_err_idx);
%             fg_convergence(i-1) = norm(fg_diff);
            fg_convergence(i-1) = sqrt(mean(fg_diff.^2,'all'));
        end
    end    
end

if ~isempty(dof_sg)
    [dof_sg,ii] = sort(dof_sg);
    error_sg = error_sg(ii);
    sg_files = sg_files(ii);
    sg_f = sg_f(ii);
    sg_convergence = [];
    if numel(sg_f)>1
        for i=2:numel(sg_f)
            sg_diff = sg_f{i}-sg_f{i-1};
%             sg_diff = sg_f{i}(z_err_idx,p_err_idx)-sg_f{i-1}(z_err_idx,p_err_idx);
%             sg_convergence(i-1) = norm(sg_diff); 
            sg_convergence(i-1) = sqrt(mean(sg_diff.^2,'all'));
        end
    end
end

figure
if ~isempty(dof_sg)
    semilogy(dof_sg(2:end),sg_convergence,'-o')
    xlabel('deg');
    ylabel('RMS(f(deg)-f(deg-1))');
    set(gca,'fontsize', 18);
    hold on
end
if ~isempty(dof_fg)
      semilogy(dof_fg(2:end),fg_convergence,'-o');
      xlabel('deg');
      ylabel('RMS(f(deg)-f(deg-1))');
      set(gca,'fontsize', 18);
end
hold off

figure
if ~isempty(dof_sg)
    loglog(dof_sg,error_sg,'-o','LineWidth',2)
    hold on
end
xlabel('Degrees of Freedom');
ylabel('RMS Absolute Error');
if ~isempty(dof_fg)
    loglog(dof_fg,error_fg,'-o','LineWidth',2)
end
hold off
set(gca,'FontSize',16)
dof_legend_strings{1} = ['Sparse Grid (deg=',num2str(opts.deg),')'];
dof_legend_strings{2} = ['Full Grid (deg=',num2str(opts.deg),')'];
legend(dof_legend_strings,'FontSize',12)


