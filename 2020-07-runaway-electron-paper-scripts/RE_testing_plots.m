root = get_root_folder();
directory = [root,'/2020-07-runaway-electron-paper-scripts/'];
file_extension = '.png';

%% Pich angle dynamics
%
% Collision term test
% -------------------
lev = 4;
deg = 1;
[err,~,~,~,~,outputs] = asgard(@fokkerplanck1_pitch_C,'timestep_method','BE','dt',0.001,...
    'num_steps',500,'lev',lev,'deg',deg,'save_output',true,'quiet',false,...
    'output_grid','dual_valued');
load(outputs.output_file_name);

figure('Position', [10 10 1600 400])
subplot(1,3,1)
plot_output1(outputs,[5],lev,deg,dt,[],0);
subplot(1,3,2)
plot_output1(outputs,[50],lev,deg,dt,[],0);
subplot(1,3,3)
plot_output1(outputs,[500],lev,deg,dt,[],0);

filename = [directory,'RE-pitch-C-deg-',num2str(deg)];
savefig([filename,'.fig']);
exportgraphics(gcf,[filename,'.pdf'],'ContentType','vector','BackgroundColor','none');

deg = 4;
[err,~,~,~,~,outputs] = asgard(@fokkerplanck1_pitch_C,'timestep_method','BE','dt',0.001,...
    'num_steps',500,'lev',lev,'deg',deg,'save_output',true,'quiet',false,...
    'output_grid','dual_valued');
load(outputs.output_file_name);

figure('Position', [10 10 1600 400])
subplot(1,3,1)
plot_output1(outputs,[5],lev,deg,dt,[],0);
subplot(1,3,2)
plot_output1(outputs,[50],lev,deg,dt,[],0);
subplot(1,3,3)
plot_output1(outputs,[500],lev,deg,dt,[],0);

filename = [directory,'RE-pitch-C-deg-',num2str(deg)];
savefig([filename,'.fig']);
exportgraphics(gcf,[filename,file_extension],'ContentType','vector','BackgroundColor','none');

% Electric field term test
% ------------------------
lev = 4;
deg = 4;
dt = 0.05;
num_steps = 50;
[err,~,~,~,~,outputs_noadapt] = asgard(@fokkerplanck1_pitch_E,...
    'timestep_method','BE','dt',dt,'num_steps',num_steps,'lev',4,'deg',4,...
    'save_output',true,'quiet',false,'adapt',false);
adapt_threshold = 1e-3;
[err_adapt,~,~,~,~,outputs_adapt] = asgard(@fokkerplanck1_pitch_E,...
    'timestep_method','BE','dt',dt,'num_steps',num_steps,'lev',4,'deg',4,...
    'save_output',true,'quiet',false,'adapt',true,...
    'adapt_threshold',adapt_threshold,'adapt_initial_condition',true,...
    'output_grid','elements');

plot_times = [10,50];
figure('Position', [10 10 1600 400])
hold off

% plot non-adapted
outputs = outputs_noadapt;
subplot(1,2,1)
plot_output1(outputs,plot_times,lev,deg,dt,adapt_threshold,1);
title('Uniform')
% plot adapted
outputs = outputs_adapt;
subplot(1,2,2)
plot_output1(outputs,plot_times,lev,deg,dt,adapt_threshold,1);
title('Adapted')

filename = [directory,'RE-pitch-E-deg-',num2str(deg)];
savefig([filename,'.fig']);
exportgraphics(gcf,[filename,file_extension],'ContentType','vector','BackgroundColor','none');

% plot the hierarchical grid
figure('Position', [10 10 1600 400]);
output = load(outputs_noadapt.output_file_name);
subplot(1,2,1)
pde = output.pde;
opts = output.opts;
hash_table = output.hash_table;
coordinates = plot_adapt(pde,opts,hash_table);
set(gca,'FontSize',14)
title('Uniform')

output = load(outputs_adapt.output_file_name);
subplot(1,2,2)
pde = output.pde;
opts = output.opts;
hash_table = output.hash_table;
coordinates = plot_adapt(pde,opts,hash_table);
set(gca,'FontSize',14)
title('Adapted')

filename = [directory,'RE-pitch-E-mesh-hierarchy-deg-',num2str(deg)];
savefig([filename,'.fig']);
exportgraphics(gcf,[filename,file_extension],'ContentType','vector','BackgroundColor','none');

% case 2 - advected pulse
lev = 4;
deg = 4;
dt = 0.04;
num_steps = 40;
adapt_threshold = 1e-2;
[err_adapt,~,~,~,~,outputs_noadapt] = asgard(@fokkerplanck1_pitch_E,'lev',lev,'deg',deg,'dt',dt,...
    'num_steps',num_steps,'case',2,'timestep_method','matrix_exponential',...
    'adapt',false,'adapt_initial_condition',true,...
    'adapt_threshold',1e-2,'output_grid','quadrature','save_output',true);

[err_adapt,~,~,~,~,outputs_adapt] = asgard(@fokkerplanck1_pitch_E,'lev',lev,'deg',deg,'dt',dt,...
    'num_steps',num_steps,'case',2,'timestep_method','matrix_exponential',...
    'adapt',true,'adapt_initial_condition',true,...
    'adapt_threshold',1e-2,'output_grid','elements','save_output',true);

plot_times = [5,10,40];
figure('Position', [10 10 1600 400])
hold off

% plot non-adapted
outputs = outputs_noadapt;
subplot(1,2,1)
plot_output1(outputs,plot_times,lev,deg,dt,adapt_threshold,0);
title('Uniform')
% plot adapted
outputs = outputs_adapt;
subplot(1,2,2)
plot_output1(outputs,plot_times,lev,deg,dt,adapt_threshold,0);
title('Adapted')

filename = [directory,'RE-pitch-E-case-2-deg-',num2str(deg)];
savefig([filename,'.fig']);
exportgraphics(gcf,[filename,file_extension],'ContentType','vector','BackgroundColor','none');

% plot the hierarchical grid
figure('Position', [10 10 1600 400]);
output = load(outputs_noadapt.output_file_name);
subplot(1,2,1)
pde = output.pde;
opts = output.opts;
hash_table = output.hash_table;
coordinates = plot_adapt(pde,opts,hash_table);
set(gca,'FontSize',14)
title('Uniform')

output = load(outputs_adapt.output_file_name);
subplot(1,2,2)
pde = output.pde;
opts = output.opts;
hash_table = output.hash_table;
coordinates = plot_adapt(pde,opts,hash_table);
set(gca,'FontSize',14)
title('Adapted')

filename = [directory,'RE-pitch-E-case-2-mesh-hierarchy-deg-',num2str(deg)];
savefig([filename,'.fig']);
exportgraphics(gcf,[filename,file_extension],'ContentType','vector','BackgroundColor','none');

disp('END');

function plot_output1(outputs,plot_times,lev,deg,dt,adapt_threshold,use_log)
green = '#77AC30';
blue  = '#0072BD';
red = '#D95319';
aqua = '#4DBEEE';
orange = '#EDB120';
purple = '#7E2F8E';
colors = {green,blue,red,aqua,orange,purple};

output = load(outputs.output_file_name);
pde = output.pde;
opts = output.opts;
hash_table = output.hash_table;
f_realspace_analytic_nD_t = output.f_realspace_analytic_nD_t;
f_realspace_nD_t = output.f_realspace_nD_t;
fval_t = output.fval_t;
nodes_t = output.nodes_t;

cnt = 1;
for i=1:numel(plot_times)
    hold on
    p=plot(nodes_t{plot_times(i)}{1},f_realspace_analytic_nD_t{plot_times(i)},'-','LineWidth',5.5,'Color',colors{i});
    p.Color(4)=0.5;
    plot(nodes_t{plot_times(i)}{1},f_realspace_nD_t{plot_times(i)},'-o','LineWidth',1,'Color',colors{i});
    labels{cnt} = ['t = ',num2str(plot_times(i)*outputs.dt),' (analytic)'];
    if opts.adapt
        labels{cnt+1} = ['t = ',num2str(plot_times(i)*outputs.dt),...
            ' (numeric), adapt threshold=',num2str(adapt_threshold),...
            ', deg=',num2str(deg),...
            ', DoF=',num2str(numel(fval_t{plot_times(i)})),...
            ', L2 err=',num2str(outputs.err{plot_times(i)},'%4.1e')];
    else
        labels{cnt+1} = ['t = ',num2str(plot_times(i)*outputs.dt),...
            ' (numeric), lev=',num2str(lev),...
            ', deg=',num2str(deg),...
            ', DoF=',num2str(numel(fval_t{plot_times(i)})),...
            ', L2 err=',num2str(outputs.err{plot_times(i)},'%4.1e')];
    end
    legend(labels,'FontSize',12,'Location','northwest');
    xlabel('\zeta');
    ylabel('$f(\zeta)$','Interpreter','latex');
    hold off
    set(gca,'FontSize',14)
    if use_log
    set(gca,'yscale','log')
    end
    cnt = cnt+2;
end
end
