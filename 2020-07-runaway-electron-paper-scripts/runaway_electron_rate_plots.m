% Produce the runaway electron fraction plots

% Calculate the relativistic and non-relativistic runaway electron rates
% from what is in Franz's thesis (pg 71 - which comes from Connor and Hastie, Nucl.
% Fusion, 15:415, 1975).

clear all;

root_folder = get_root_folder();

set(groot,'defaultLineLineWidth',2.0);
set(groot,'defaultAxesFontSize',20);
set(groot,'defaultLegendFontSize',20);

red = '#E74C3C';
blue = '#2E86C1';
green = '#229954';

cgs.e = 4.803e-10;
cgs.m_e = 9.109e-28;
cgs.c = 2.997e10;
cgs.amu = 1.66e-24;
cgs.kB = 1.38e-16;

SI.e = 1.602e-19;
SI.c = 3e8;
SI.m_e = 9.1e-31;
SI.amu = 1.66e-27;
SI.eps0 = 8.854e-12;
SI.kB = 1.38e-23;

Z_b = 1;
energy_eV = 20e3;
T_e_eV = 2/3*energy_eV;
T_e_K = T_e_eV * 11605;

% calculate Dreicer E field in cgs

cgs.n_e = 5e13;
cgs.coulomb_ln = 23.5-log(cgs.n_e^(1/2)*T_e_K^(-5/4))-sqrt(10^-5+((log(T_e_K)-2)^2)/16); % NRL forumulary
coulomb_ln = cgs.coulomb_ln;
cgs.th_e = cgs.m_e * cgs.c^2 / (T_e_eV);
% cgs.v_th = 4.19e7*T_e_eV^(1/2);
cgs.v_th = sqrt(cgs.kB*T_e_K/cgs.m_e);
cgs.E_D = 4*pi * cgs.e^3 * cgs.n_e * cgs.coulomb_ln / (cgs.m_e * cgs.v_th^2); % from Connor Hastie 1975
cgs.E_R = 4*pi*cgs.n_e*cgs.e^3*cgs.coulomb_ln/(cgs.m_e*cgs.c^2);
cgs.v_c = @(E_cgs) sqrt(4*pi*cgs.e^3 * cgs.n_e * cgs.coulomb_ln / (cgs.m_e * E_cgs));
disp(['E_D V/m (from cgs): ',num2str(cgs.E_D*300*100)]);
disp(['E_R V/m (from cgs): ',num2str(cgs.E_R*300*100)]);
disp(['E_R V/m (from cgs 2): ',num2str(cgs.E_D*300*100*cgs.kB*T_e_K/(cgs.m_e*cgs.c^2))]);
disp(['E_D / E_R : ', num2str(cgs.E_D/cgs.E_R)]);
disp(['v_th m/s (from cgs): ',num2str(cgs.v_th/100)]);

% calculate Dreicer E field in SI

SI.n_e = cgs.n_e*1e6;
SI.v_th = sqrt(T_e_eV*SI.e/SI.m_e);
SI.coulomb_ln = coulomb_ln;
SI.E_D = 1/(4*pi*SI.eps0^2) * SI.n_e * SI.e^3 * SI.coulomb_ln / (SI.m_e * SI.v_th^2);
disp(['E_D V/m (from SI): ', num2str(SI.E_D)]);
disp(['v_th m/s (from SI): ',num2str(SI.v_th)]);

% calculate the RE production rate

cgs.nu_ee2 = 4*pi * cgs.n_e * cgs.e^4 * cgs.coulomb_ln / (cgs.th_e^(3/2) * cgs.m_e^3 * cgs.c^3);
cgs.nu_ee = 4*pi * cgs.n_e * cgs.e^4 * cgs.coulomb_ln / (cgs.m_e^2 * cgs.v_th^3);
disp(['nu_ee 1/s (from cgs): ',num2str(cgs.nu_ee)]);
disp(['nu_ee 1/s (from cgs 2): ',num2str(cgs.nu_ee2)]);

% normalized (norm) parameters for use in FokkerPlanck2_complete_div
% here we have 3 types of variables ...
% 1) T,n,v_th - absolute values in SI
% 2) T_ref, n_ref, v_th_ref - these are the "reference values" (or the "tilde" values in the manuscript)
% 3) T_norm, n_norm, v_th_norm - these are T_norm = T / T_ref

ref.energy_eV = 20e3;
ref.T_eV = 2/3*ref.energy_eV;
ref.T_K = ref.T_eV * 11605;
ref.n_e_m3 = 5e19;
ref.v_th = sqrt (2*ref.T_eV*SI.e / SI.m_e);
ref.coulomb_ln = 23.5-log((ref.n_e_m3*1e-6)^(1/2)*ref.T_K^(-5/4))-sqrt(10^-5+((log(ref.T_K)-2)^2)/16);
ref.E_D = SI.e^3 * ref.n_e_m3 * ref.coulomb_ln / (4*pi*SI.eps0^2 * ref.T_eV * SI.e);

norm.delta = SI.v_th / SI.c;
norm.v_th = sqrt(T_e_eV / ref.T_eV);
norm.p = SI.m_e * ref.v_th;
norm.n = SI.n_e / ref.n_e_m3; %SI.n_e;
norm.nu_ee = (ref.T_eV/T_e_eV)^(3/2) * coulomb_ln / ref.coulomb_ln;
norm.tau = 10^5; % this is relative to the B field or some other thing not present in these calcs
norm.E = ref.E_D/2;
norm.n = SI.n_e/ref.n_e_m3;

disp(['delta: ', num2str(norm.delta)]);
disp(['tau: ', num2str(norm.tau)]); % this is relative to the B field or some other thing not present in these calcs
disp(['E: ', num2str(norm.E), '  (this is E_D/2)']);
disp(['v_th: ', num2str(norm.v_th)]);
disp(['nu_ee: ', num2str(norm.nu_ee)]);
disp(['n: ', num2str(norm.n)]);

% Scan over E/E_D (0 to 0.6) for Zb=1 - compare with Fig 9 (pg 68) of Franz

ratio_max = 0.2;
N = 100;
ratio = linspace(0.,ratio_max,N);

% Franz
Cn = 1/SI.n_e;
T1 = @(E,Z_b) Cn .* cgs.nu_ee .* (E ./ cgs.E_D).^(3/16*Z_b+1);
T2 = @(E,Z_b) exp(-cgs.E_D./(4.*E)-sqrt((1+Z_b).*cgs.E_D./E));

% Kruskal & Bernstein
kruskal_bernstein = @(E) 0.35 .* E.^(-3/8) .* exp(-((2./E).^(1/2)+1./(4*E)));

% Connor & Hastie 1975 eq (67)
C = 1;
T1 = @(E,Z) C .* cgs.n_e * cgs.nu_ee * (E/cgs.E_D).^(-3/16*(Z+1));
T2 = @(E,Z) exp(-cgs.E_D./(4*E)-sqrt((1+Z).*cgs.E_D./E));
T3 = @(E,Z) exp( -cgs.kB*T_e_K./(cgs.m_e.*cgs.c^2) * ( 1/8*(cgs.E_D./E).^2 + 2/3*(cgs.E_D./E).^(3/2).*(1+Z).^(1/2) )  );
% T3 = @(E,Z) exp( -norm.delta.^2 * 0.5 * ( 1/8*(cgs.E_D./E).^2 + 2/3*(cgs.E_D./E).^(3/2).*(1+Z).^(1/2) )  ); % account for the vth^2=2T/m vs vth^2=T/m difference

connor_hastie_nr = @(E,Z) T1(E,Z) .* T2(E,Z);
connor_hastie_r =  @(E,Z) connor_hastie_nr(E,Z) .* T3(E,Z);

% Kulsrud data

kulsrud_E1 = [0.04,0.06,0.08,0.1];
kulsrud_Z1 = [1.914e-6,5.411e-5,3.177e-4,1.004e-3];

kulsrud_E2 = [0.06,0.08,0.1];
kulsrud_Z2 = [2.611e-5,1.735e-4,5.839e-4];

kulsrud_E3 = [0.08,0.1];
kulsrud_Z3 = [1.047e-4,3.757e-4];

kulsrud_E10 = [0.08,0.1];
kulsrud_Z10 = [9.0e-6,4.49e-5];

kulsrud_E008 = [1,2,3,10];
kulsrud_Z008 = [3.177e-4,1.735e-4,1.047e-4,9.0e-6];

CH_normfac_nr = 1/(connor_hastie_nr(0.1*cgs.E_D,1)/kulsrud_Z1(end));
CH_normfac_r = CH_normfac_nr;
KB_normfac = 1/(kruskal_bernstein(0.1)/kulsrud_Z1(end));

args.p_max = 10;
args.v_th = norm.v_th;
args.nu_ee = norm.nu_ee;
args.tau = norm.tau;
args.n = norm.n;

% Scan over E/E_D ratio
do_E_scan = true;
if do_E_scan
    
    N2 = 15;
    ratio_max2 = 0.1;
    ratio2 = linspace(0.02,ratio_max2,N2);
    num_steps = 10;
    deg = 4;
    lev = 4;
    
    alphas_filename = [root_folder,'/output/alpha_E_scan.mat'];
    load_alphas = true;
    if load_alphas
        load(alphas_filename);
    else
        
        % generate or load asgard run output
        
        get_file_id = @(str,E,Z) [str,num2str(E,'%1.3f'),'_Z',num2str(Z)];
        
        % non-relativistic
        args.delta = 0;
        for Z = [1 2 10]
            for i=1:N2
                args.E = ratio2(i) * SI.E_D / norm.E; % E into asgard is 2*E/E_D per the normalization
                args.Z = Z;
                dt = 2./args.E.^2/num_steps;
                if dt > 2000
                    dt = 2000;
                end
                file_id = get_file_id('E_scan_nonrel_E',args.E,args.Z);
                [~,~,~,~,~,~,opts] = asgard(@fokkerplanck2_complete,...
                    'timestep_method','BE','num_steps',num_steps,'dt',dt,'deg',deg,'lev',lev,'case',5,...
                    'cmd_args',args,'quiet',true,'calculate_mass',true,'grid_type','SG',...
                    'update_params_each_timestep',true,'save_output',true,'output_filename_id',file_id,...
                    'save_freq',num_steps);
                output_filename = create_output_filename(opts);
                load(output_filename,'outputs');
                alpha_nr(i,Z) = outputs.alpha_t{end};
                clear outputs
            end
        end
        
        % relativistic
        args.delta = norm.delta;
        for Z = [1 2 10]
            for i=1:N2
                args.E = ratio2(i) * SI.E_D / norm.E; % E into asgard is 2*E/E_D per the normalization
                args.Z = Z;
                dt = 2./args.E.^2/num_steps;
                if dt > 2000
                    dt = 2000;
                end
                file_id = get_file_id('E_scan_relati_E',args.E,args.Z);
                [~,~,~,~,~,~,opts] = asgard(@fokkerplanck2_complete,...
                    'timestep_method','BE','num_steps',num_steps,'dt',dt,'deg',deg,'lev',lev,'case',5,...
                    'cmd_args',args,'quiet',true,'calculate_mass',true,'grid_type','SG',...
                    'update_params_each_timestep',true,'save_output',true,'output_filename_id',file_id,...
                    'save_freq',num_steps);
                output_filename = create_output_filename(opts);
                load(output_filename,'outputs');
                alpha_r(i,Z) = outputs.alpha_t{end};
                clear outputs
            end
        end
        
        save(alphas_filename,'alpha_r','alpha_nr');
        
    end
    
    % plot the results for the loaded data
    
    set(groot,'defaultLineLineWidth',2.0)
    set(groot,'defaultAxesFontSize',20)
    set(groot,'defaultLegendFontSize',14)
    figure('Position',[0 0 600 600])
    semilogy(ratio,kruskal_bernstein(ratio)*KB_normfac,'LineWidth',10,'Color','#E8E7E7','DisplayName','K-B (nr)');
    title('RE production rate benchmark','Interpreter','latex');
    xlabel('$E/E_D$','Interpreter','latex');
    ylabel('RE production rate $\alpha$ [arb. units]','Interpreter','latex');
    ylim([1e-10 1]);
    hold on
    
    semilogy(ratio,connor_hastie_nr(ratio*cgs.E_D,1)*CH_normfac_nr,'LineStyle','--','DisplayName','C-H (nr, Z=1)','color',red);
    semilogy(ratio,connor_hastie_r(ratio*cgs.E_D,1)*CH_normfac_r,'LineStyle',':','DisplayName','C-H (r, Z=1)','color',red);
    semilogy(ratio,connor_hastie_nr(ratio*cgs.E_D,10)*CH_normfac_nr,'LineStyle','--','DisplayName','C-H (nr, Z=10)','color',blue);
    semilogy(ratio,connor_hastie_r(ratio*cgs.E_D,10)*CH_normfac_r,'LineStyle',':','DisplayName','C-H (r, Z=10)','color',blue);
    semilogy(kulsrud_E1,kulsrud_Z1,'Marker','s','MarkerSize',26,'MarkerFaceColor','auto','DisplayName','Kulsrud (nr, Z=1)','color',red);
    semilogy(kulsrud_E2,kulsrud_Z2,'Marker','s','MarkerSize',26,'MarkerFaceColor','auto','DisplayName','Kulsrud (nr, Z=2)','color',green);
    semilogy(kulsrud_E10,kulsrud_Z10,'Marker','s','MarkerSize',26,'MarkerFaceColor','auto','DisplayName','Kulsrud (nr, Z=10)','color',blue);
    norm_fac = kulsrud_Z1(end)/alpha_nr(N2,1);
    disp(['norm_fac: ',num2str(norm_fac)]);
    hold on
    semilogy(ratio2,alpha_nr(:,1)*norm_fac,'LineWidth',5,'DisplayName','ASGarD (nr, Z=1)','color',red);
    semilogy(ratio2,alpha_nr(:,2)*norm_fac,'LineWidth',5,'DisplayName','ASGarD (nr, Z=2)','color',green);
    semilogy(ratio2,alpha_nr(:,10)*norm_fac,'LineWidth',5,'DisplayName','ASGarD (nr, Z=10)','color',blue);
    
    semilogy(ratio2,alpha_r(:,1)*norm_fac,'LineWidth',5,'DisplayName','ASGarD (r, Z=1)','color',red,'LineStyle',':');
    semilogy(ratio2,alpha_r(:,2)*norm_fac,'LineWidth',5,'DisplayName','ASGarD (r, Z=2)','color',green,'LineStyle',':');
    semilogy(ratio2,alpha_r(:,10)*norm_fac,'LineWidth',5,'DisplayName','ASGarD (r, Z=10)','color',blue,'LineStyle',':');
    legend('FontSize',10)
    hold off
end

% Compare DoF required for given accuracy in alpha

args.delta = 0;
args.Z = 1;
deg_span_FG = [2,3,4];
lev_span_FG = [4,5];
num_steps = 20;
E = [0.1,0.07,0.055];
E_string = {'010','007','0055'};


do_dof_scan = true;
if do_dof_scan
    
    get_file_id = @(str,Estr,lev,deg) [str,Estr,'_lev',num2str(lev),'_deg',num2str(deg)];
    
    % generate (or load) the FG reference data (both to demonstrate a convergence in
    % the FG results, and to use as a "correct" answer to compare SG and
    % ASG to).
    
    alphas_filename = [root_folder,'/output/alpha_FG.mat'];
    load_alphas = true;
    if load_alphas
        load(alphas_filename);
    else
        load(alphas_filename); % load so we can add to it
        for i=1:numel(E)
            args.E = E(i);
            dt = 2./args.E.^2/num_steps;
            
            for lev=lev_span_FG
                for deg=deg_span_FG
                    if isempty(alpha_FG{i,lev,deg})
                        file_id = get_file_id('FG_',E_string{i},lev,deg)
                        [~,~,~,~,~,~,opts] = asgard(@fokkerplanck2_complete,...
                            'timestep_method','BE','num_steps',num_steps,'dt',dt,'deg',deg,'lev',lev,'case',5,...
                            'cmd_args',args,'quiet',false,'calculate_mass',true,'grid_type','FG',...
                            'update_params_each_timestep',true,'output_grid','fixed',...
                            'save_output',true,'output_filename_id',file_id,'time_independent_build_A',false,...
                            'save_freq',num_steps);
                        output_filename = create_output_filename(opts);
                        load(output_filename,'outputs');
                        alpha_FG{i,lev,deg} = outputs.alpha_t;
                        dof_FG{i,lev,deg} = numel(outputs.fval_t{end});
                        time_FG{i,lev,deg} = outputs.time_array;
                        clear outputs
                        save(alphas_filename,'alpha_FG','dof_FG','time_FG');
                    end
                end
            end
            
            
        end
        
    end
    
    
    % plot the convergence of the FG results
    
    set(groot,'defaultLineLineWidth',2.0)
    set(groot,'defaultAxesFontSize',20)
    set(groot,'defaultLegendFontSize',14)
    ms = 12;
    lw = 2;
    linestyles = {'--',':','-.','-','--'};
    linecolors = {red,blue,green};
    
    plot_FG_convergence = true;
    if plot_FG_convergence
        
        figure('Position',[0 0 600 600])
        
        for i=1:numel(E)
            args.E = E(i);
            
            for deg=4:7
                alpha_t = alpha_FG{i,4,deg};
                dof = dof_FG{i,4,deg};
                this_lw = lw;
                if (i==1&&deg==4) || (i==2&&deg==5) || (i==3&&deg==6)
                    this_lw=6;
                end
                p = semilogy(time_FG{i,4,deg}(2:end),cell2mat(alpha_FG{i,4,deg}),'DisplayName',['FG (E/E_D=',num2str(E(i),'%1.3f'),' deg=',num2str(deg),',lev=',num2str(lev), ', DoF=',num2str(dof)],'LineWidth',this_lw,'Color',linecolors{i},'LineStyle',linestyles{deg-4+1});
                hold off
                if (i==1&&deg==4) || (i==2&&deg==5) || (i==3&&deg==6)
                    p.Color(4) = 0.5;
                end
                hold on
            end
            
        end
        xlabel('{Normalized Time ($\tilde{t}=t/\tilde{\nu}_{ee}$)}','Interpreter','latex');
        ylabel('{RE Production Rate ($\alpha(\tilde{t}))$ [arb. units]}','Interpreter','latex');
        title('Full Grid (FG) Convergence of $\alpha$','Interpreter','latex');
        legend('FontSize',10)
        ylim([10^-8 10^-2])
        xlim([0 500])
        legend
        hold off
        
    end
    
    
    % now look at how the SG solution converges to the "correct" solution
    
    deg_span_SG = [2,3,4];
    lev_span_SG = [4,5,6];
    
    alphas_filename = [root_folder,'/output/alpha_SG.mat'];
    load_alphas = true;
    if load_alphas
        load(alphas_filename);
    else
        
        for i=1:numel(E)
            args.E = E(i);
            dt = 2./args.E.^2/num_steps;
            
            for lev=lev_span_SG
                for deg=deg_span_SG
                    file_id = get_file_id('SG_',E_string{i},lev,deg);
                    [~,~,~,~,~,~,opts] = asgard(@fokkerplanck2_complete,...
                        'timestep_method','BE','num_steps',num_steps,'dt',dt,'deg',deg,'lev',lev,'case',5,...
                        'cmd_args',args,'quiet',true,'calculate_mass',true,'grid_type','SG',...
                        'update_params_each_timestep',true,'output_grid','fixed',...
                        'save_output',true,'output_filename_id',file_id,'time_independent_build_A',true,...
                        'save_freq',num_steps);
                    output_filename = create_output_filename(opts);
                    load(output_filename,'outputs');
                    alpha_SG{i,lev,deg} = outputs.alpha_t;
                    dof_SG{i,lev,deg} = numel(outputs.fval_t{end});
                    time_SG{i,lev,deg} = outputs.time_array;
                    clear outputs
                end
            end
            
        end
        
        save(alphas_filename,'alpha_SG','dof_SG','time_SG');
        
    end
    
    plot_SG_convergence = true;
    if plot_SG_convergence
        
        load([root_folder,'/output/alpha_FG.mat']);
        
        for deg=deg_span_SG
            
            for lev=4:5
                
                figure('Position',[0 0 600 600])
                
                for i=1:numel(E)
                    args.E = E(i);
                    
                    % include "correct" FG results for comparison
                    
                    this_lw = lw;
                    switch i
                        case 1
                            deg_FG=4;
                        case 2
                            deg_FG=5;
                        case 3
                            deg_FG=6;
                    end
                    alpha_t = alpha_FG{i,4,deg_FG};
                    dof = dof_FG{i,4,deg_FG};
                    p = semilogy(time_FG{i,4,deg_FG}(2:end),cell2mat(alpha_FG{i,4,deg_FG}),'DisplayName',['Converged Solution (E/E_D=',num2str(E(i),'%1.3f'),' FG deg=',num2str(deg_FG),',lev=',num2str(lev), ', DoF=',num2str(dof)],'LineWidth',6,'Color',linecolors{i},'LineStyle','-');
                    hold on
                    p.Color(4) = 0.5;
                    
                    % now overlay SG results
                    alpha_t = alpha_SG{i,lev,deg};
                    dof = dof_SG{i,lev,deg};
                    this_lw = lw;
                    this_alpha = cell2mat(alpha_SG{i,lev,deg});
                    this_alpha(this_alpha<1e-8)=1e-8; % just for plotting purposes on the log scale
                    p = semilogy(time_SG{i,lev,deg}(2:end),this_alpha,'DisplayName',['SG (E/E_D=',num2str(E(i),'%1.3f'),' deg=',num2str(deg),',lev=',num2str(lev), ', DoF=',num2str(dof)],'LineWidth',this_lw,'Color',linecolors{i},'LineStyle','--');
                    if lev==4 || lev==5% overlay FG
                        alpha_t = alpha_FG{i,lev,deg};
                        dof = dof_FG{i,lev,deg};
                        p = semilogy(time_FG{i,lev,deg}(2:end),cell2mat(alpha_FG{i,lev,deg}),'DisplayName',['FG (E/E_D=',num2str(E(i),'%1.3f'),' deg=',num2str(deg),',lev=',num2str(lev), ', DoF=',num2str(dof)],'LineWidth',this_lw,'Color',linecolors{i},'LineStyle','-');
                    end
                    
                end
                
                xlabel('{Normalized Time ($\tilde{t}=t/\tilde{\nu}_{ee}$)}','Interpreter','latex');
                ylabel('{RE Production Rate ($\alpha(\tilde{t}))$ [arb. units]}','Interpreter','latex');
                title(['Sparse Grid (SG) Convergence of $\alpha$ for deg=',num2str(deg),' lev=',num2str(lev)],'Interpreter','latex');
                legend('FontSize',10)
                ylim([10^-8 10^-2])
                xlim([0 500])
                legend
                hold off
                
            end
            
        end
        
    end
    
    % now look at how the ASG solution converges to the "correct" solution
    
    deg_span_ASG = [3];
    lev_span_ASG = [5,6]; % lev for the ASG refers to the adaptivity threshold, i.e., 10^-lev
    lev = 3;
    
    alphas_filename = [root_folder,'/output/alpha_ASG.mat'];
    load_alphas = false;
    if load_alphas
        load(alphas_filename);
    else
        if exist(alphas_filename,'file')            
            load(alphas_filename);
        end
        for i=1:numel(E)
            args.E = E(i);
            dt = 2./args.E.^2/num_steps;
                        
            for lev=lev_span_ASG
                for deg=deg_span_ASG
                    empty = true;
                    try
                        if ~isempty(alpha_ASG{i,lev,deg});empty=false;end
                    catch
                    end
                    if empty
                        file_id = get_file_id('ASG_',E_string{i},lev,deg);
                        [~,~,~,~,~,~,opts] = asgard(@fokkerplanck2_complete,'timestep_method','BE','num_steps',num_steps,...
                            'dt',dt,'deg',deg,'lev',3,'case',5,'cmd_args',args,'quiet',false,'calculate_mass',true,...
                            'grid_type','SG','update_params_each_timestep',true,'output_grid','fixed','save_output',true,...
                            'output_filename_id',file_id,'time_independent_build_A',false,...
                            'adapt',true,'adapt_initial_condition',true,'adapt_threshold',1*10^(-lev),...
                            'max_lev',8,'max_lev_coeffs',true,'save_freq',num_steps);
                        output_filename = create_output_filename(opts);
                        load(output_filename,'outputs');
                        alpha_ASG{i,lev,deg} = outputs.alpha_t;
                        dof_ASG{i,lev,deg} = numel(outputs.fval_t{end});
                        time_ASG{i,lev,deg} = outputs.time_array;
                        clear outputs
                        save(alphas_filename,'alpha_ASG','dof_ASG','time_ASG');
                    end
                end
            end
            
        end
    end
    
    plot_ASG_convergence = true;
    if plot_ASG_convergence
        
        load([root_folder,'/output/alpha_FG.mat']);
        
        for deg=deg_span_ASG
            
            
            figure('Position',[0 0 600 600])
            
            for i=1:numel(E)
                args.E = E(i);
                
                % include "correct" FG results for comparison
                
                this_lw = lw;
                switch i
                    case 1
                        deg_FG=4;
                    case 2
                        deg_FG=5;
                    case 3
                        deg_FG=6;
                end
                alpha_t = alpha_FG{i,4,deg_FG};
                dof = dof_FG{i,4,deg_FG};
                p = semilogy(time_FG{i,4,deg_FG}(2:end),cell2mat(alpha_FG{i,4,deg_FG}),'DisplayName',['Converged Solution (E/E_D=',num2str(E(i),'%1.3f'),' FG deg=',num2str(deg_FG),',lev=',num2str(lev), ', DoF=',num2str(dof)],'LineWidth',6,'Color',linecolors{i},'LineStyle','-');
                hold on
                p.Color(4) = 0.5;
                
                % now overlay ASG results
                for lev=lev_span_ASG
                    alpha_t = alpha_ASG{i,lev,deg};
                    dof = dof_ASG{i,lev,deg};
                    this_lw = lw;
                    this_alpha = cell2mat(alpha_ASG{i,lev,deg});
                    this_alpha(this_alpha<1e-8)=1e-8; % just for plotting purposes on the log scale
                    p = semilogy(time_ASG{i,lev,deg}(2:end),this_alpha,...
                        'DisplayName',['ASG (E/E_D=',num2str(E(i),'%1.3f'),' deg=',num2str(deg),...
                        ', ASG Threshold=',num2str(1*10^(-lev),'%1.0e'), ...
                        ', DoF=',num2str(dof)],'LineWidth',this_lw,'Color',linecolors{i},'LineStyle',linestyles{lev-2+1});
                end
                
            end
            
            xlabel('{Normalized Time ($\tilde{t}=t/\tilde{\nu}_{ee}$)}','Interpreter','latex');
            ylabel('{RE Production Rate ($\alpha(\tilde{t}))$ [arb. units]}','Interpreter','latex');
            title(['{Adaptive Sparse Grid (ASG) Convergence of $\alpha$ for deg=',num2str(deg),'}'],'Interpreter','latex');
            legend('FontSize',10)
            ylim([10^-8 10^-2])
            xlim([0 500])
            legend
            hold off
            
        end
        
    end
        
end
disp('');

