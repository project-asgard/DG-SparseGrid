root = get_root_folder();
directory = [root,'/2020-07-runaway-electron-paper-scripts/'];
file_extension = '.png';


%% Momentum dynamics

% Momentum (Collision, E=R=0)
% Two initial conditions:
% Case 1: smooth Maxwellian
% Case 2: discontinuous step function

lev = 5;
deg = 4; 
dt = 50;
[err,~,~,~,~,outputs] = asgard(@fokkerplanck1_momentum_C_div,...
    'timestep_method','matrix_exponential','lev',lev,'deg',deg,'num_steps',5,'dt',dt,...
    'case',2,'save_output',true,'calculate_mass',true,'quiet',true);
load(outputs.output_file_name);
figure('Position', [10 10 1600 400])
subplot(1,3,1)
plot_times = [];
plot_output1(outputs,plot_times,lev,deg,dt,[],1);

% look at the error versus lev
dt = 350;
for deg=1:4
for lev=3:6
    disp(['Solving for deg=',num2str(deg)]);
    disp(['Solving for lev=',num2str(lev)]);
  [err,~,~,~,~,outputs] = asgard(@fokkerplanck1_momentum_C_div,...
    'timestep_method','matrix_exponential','lev',lev,'deg',deg,'num_steps',1,'dt',dt,...
    'case',2,'save_output',true,'calculate_mass',true,'quiet',true);  
    err_table(lev,deg) = err;
end
end

%% Pitch angle dynamics

% Pitch (Collision, E=R=0)

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

function plot_output1(outputs,plot_indicies,lev,deg,dt,adapt_threshold,use_log)

if isempty(plot_indicies)
    num_times = numel(outputs.fval_t);
    plot_indicies = linspace(1,num_times,num_times);
end

green = '#77AC30';
blue  = '#0072BD';
red = '#D95319';
aqua = '#4DBEEE';
orange = '#EDB120';
purple = '#7E2F8E';
colors = {green,blue,red,aqua,orange,purple};

loaded_struct = load(outputs.output_file_name);
pde = loaded_struct.pde;
opts = loaded_struct.opts;
hash_table = loaded_struct.hash_table;
f_realspace_analytic_nD_t = loaded_struct.outputs.f_realspace_analytic_nD_t;
f_realspace_nD_t = loaded_struct.outputs.f_realspace_nD_t;
fval_t = loaded_struct.outputs.fval_t;
nodes_t = loaded_struct.outputs.nodes_t;

cnt = 1;
for i=1:numel(plot_indicies)
    hold on
    color_index = 1+mod(i-1,numel(colors));
    p=plot(nodes_t{plot_indicies(i)}{1},f_realspace_analytic_nD_t{plot_indicies(i)},'-','LineWidth',5.5,'Color',colors{color_index});
    p.Color(4)=0.5;
    plot(nodes_t{plot_indicies(i)}{1},f_realspace_nD_t{plot_indicies(i)},'-o','LineWidth',1,'Color',colors{color_index});
    labels{cnt} = ['t = ',num2str((plot_indicies(i)-1)*outputs.dt),' (analytic)'];
    if opts.adapt
        labels{cnt+1} = ['t = ',num2str(plot_indicies(i)*outputs.dt),...
            ' (numeric), adapt threshold=',num2str(adapt_threshold),...
            ', deg=',num2str(deg),...
            ', DoF=',num2str(numel(fval_t{plot_indicies(i)})),...
            ', L2 err=',num2str(outputs.err{plot_indicies(i)},'%4.1e')];
    else
        labels{cnt+1} = ['t = ',num2str(plot_indicies(i)*outputs.dt),...
            ' (numeric), lev=',num2str(lev),...
            ', deg=',num2str(deg),...
            ', DoF=',num2str(numel(fval_t{plot_indicies(i)})),...
            ', L2 err=',num2str(outputs.err{plot_indicies(i)},'%4.1e')];
    end
    legend(labels,'FontSize',12,'Location','northwest');
    xlabel('\zeta');
    ylabel('$f(\zeta)$','Interpreter','latex');
    hold off
    set(gca,'FontSize',14)
    if use_log
    set(gca,'yscale','log')
    ylim([1e-10,10]);
    end
    cnt = cnt+2;
end
end
