% Produce the runaway electron fraction plots

% Calculate the relativistic and non-relativistic runaway electron rates
% from what is in Franz's thesis (pg 71 - which comes from Connor and Hastie, Nucl.
% Fusion, 15:415, 1975).

clear all;

set(groot,'defaultLineLineWidth',2.0)
set(groot,'defaultAxesFontSize',20)
set(groot,'defaultLegendFontSize',20)


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


% Scan over E/E_D ratio
do_E_scan = true;
if do_E_scan
    figure
    semilogy(ratio,kruskal_bernstein(ratio)*KB_normfac,'LineWidth',10,'Color','#E8E7E7','DisplayName','K-B (nr)');
    ylim([1e-10 1]);
    hold on
    semilogy(ratio,connor_hastie_nr(ratio*cgs.E_D,1)*CH_normfac_nr,'LineStyle','--','DisplayName','C-H (nr, Z=1)','color','red');
    semilogy(ratio,connor_hastie_r(ratio*cgs.E_D,1)*CH_normfac_r,'LineStyle',':','DisplayName','C-H (r, Z=1)','color','red');
    semilogy(ratio,connor_hastie_nr(ratio*cgs.E_D,10)*CH_normfac_nr,'LineStyle','--','DisplayName','C-H (nr, Z=10)','color','blue');
    semilogy(ratio,connor_hastie_r(ratio*cgs.E_D,10)*CH_normfac_r,'LineStyle',':','DisplayName','C-H (r, Z=10)','color','blue');
    semilogy(kulsrud_E1,kulsrud_Z1,'Marker','s','MarkerSize',26,'MarkerFaceColor','auto','DisplayName','Kulsrud (nr, Z=1)','color','red');
    semilogy(kulsrud_E2,kulsrud_Z2,'Marker','s','MarkerSize',26,'MarkerFaceColor','auto','DisplayName','Kulsrud (nr, Z=2)','color','green');
    semilogy(kulsrud_E10,kulsrud_Z10,'Marker','s','MarkerSize',26,'MarkerFaceColor','auto','DisplayName','Kulsrud (nr, Z=10)','color','blue');
    %     legend('K-B','C-H (nr)','C-H (r)','Kulsrud (Z=1)','Kulsrud (Z=2)','Kulsrud (Z=10)','ASGarD');
    legend
    
    N2 = 15;
    ratio_max2 = 0.1;
    ratio2 = linspace(0.02,ratio_max2,N2);
    args.p_max = 10;
    args.v_th = norm.v_th;
    args.nu_ee = norm.nu_ee;
    args.tau = norm.tau;
    args.n = norm.n;
    num_steps = 10;
    
    % Non-relativistic
    args.delta = 0;
    for Z = [1 2 10]
        for i=1:N2
            args.E = ratio2(i) * SI.E_D / norm.E; % E into asgard is 2*E/E_D per the normalization
            args.Z = Z;
            dt = 2./args.E.^2/num_steps;
            if dt > 2000
                dt = 2000;
            end
            [~,~,~,~,~,outputs(i)] = asgard(@fokkerplanck2_complete_div,'timestep_method','BE','num_steps',num_steps,'dt',dt,'deg',4,'lev',4,'case',5,'cmd_args',args,'quiet',true,'calculate_mass',true,'grid_type','SG','update_params_each_timestep',true);
            alpha_nr(i,Z) = outputs(i).alpha_t{end};
        end
    end
    
    % Relativistic
    args.delta = norm.delta;
    for Z = [1 2 10]
        for i=1:N2
            args.E = ratio2(i) * SI.E_D / norm.E; % E into asgard is 2*E/E_D per the normalization
            args.Z = Z;
            dt = 2./args.E.^2/num_steps;
            if dt > 2000
                dt = 2000;
            end
            [~,~,~,~,~,outputs(i)] = asgard(@fokkerplanck2_complete_div,'timestep_method','BE','num_steps',num_steps,'dt',dt,'deg',4,'lev',4,'case',5,'cmd_args',args,'quiet',true,'calculate_mass',true,'grid_type','SG','update_params_each_timestep',true);
            alpha_r(i,Z) = outputs(i).alpha_t{end};
        end
    end
    
    norm_fac = kulsrud_Z1(end)/alpha_nr(N2,1);    
    disp(['norm_fac: ',num2str(norm_fac)]);
    hold on
    semilogy(ratio2,alpha_nr(:,1)*norm_fac,'LineWidth',5,'DisplayName','ASGarD (nr, Z=1)','color','red');
    semilogy(ratio2,alpha_nr(:,2)*norm_fac,'LineWidth',5,'DisplayName','ASGarD (nr, Z=2)','color','green');
    semilogy(ratio2,alpha_nr(:,10)*norm_fac,'LineWidth',5,'DisplayName','ASGarD (nr, Z=10)','color','blue');
    
    semilogy(ratio2,alpha_r(:,1)*norm_fac,'LineWidth',5,'DisplayName','ASGarD (r, Z=1)','color','red','LineStyle',':');
    semilogy(ratio2,alpha_r(:,2)*norm_fac,'LineWidth',5,'DisplayName','ASGarD (r, Z=2)','color','green','LineStyle',':');
    semilogy(ratio2,alpha_r(:,10)*norm_fac,'LineWidth',5,'DisplayName','ASGarD (r, Z=10)','color','blue','LineStyle',':');
    legend
    hold off
end


% Scan over Z
do_Z_scan = false;
if do_Z_scan
    Z = linspace(1,10,100);
    E_ED = 0.08;
    figure
    semilogy(Z,connor_hastie_nr(E_ED*cgs.E_D,Z));
    hold on
    semilogy(Z,connor_hastie_r(E_ED*cgs.E_D,Z));
    norm_fac = connor_hastie_nr(E_ED*cgs.E_D,1)/kulsrud_Z008(1);
    semilogy(kulsrud_E008,kulsrud_Z008*norm_fac,'Marker','s','MarkerSize',16,'MarkerFaceColor','auto');
    % ylim([1e-10 1]);
    
    Z2 = 1:10;
    for i=1:numel(Z2)
        args.E = E_ED * SI.E_D / norm.E; % E into asgard is 2*E/E_D per the normalization
        args.Z = Z2(i);
        args.p_max = 10;
        args.v_th = norm.v_th;
        args.delta = 0;
        args.nu_ee = norm.nu_ee;
        args.tau = norm.tau;
        args.n = norm.n;
        
        E_cgs = args.E / 300 / 100;
        v_c_cgs = cgs.v_c(E_cgs+1e-6);
        v_c_SI = v_c_cgs/100;
        v_c_norm = v_c_SI * SI.m_e / norm.p;
        disp(['v_c (cm/s): ',num2str(v_c_cgs)]);
        disp(['v_c (m/s): ',num2str(v_c_SI)]);
        disp(['v_c (norm): ', num2str(v_c_norm)]);
        
        num_steps = 10;
        dt = 2./args.E.^2/num_steps;
        disp(i);
        [~,~,~,~,~,outputs(i)] = asgard(@fokkerplanck2_complete_div,'timestep_method','BE','num_steps',num_steps,'dt',dt,'deg',4,'lev',4,'case',5,'cmd_args',args,'quiet',true,'calculate_mass',true,'grid_type','SG','update_params_each_timestep',true);
        alpha_Z(i) = outputs(i).alpha_t{end};
    end
    
    norm_fac = connor_hastie_nr(E_ED*cgs.E_D,1)/alpha_Z(1);
    semilogy(Z2,alpha_Z*norm_fac,'Marker','none','MarkerSize',16,'MarkerFaceColor','auto');
    legend('C-H (nr)','C-H (r)','kulsrud (E/E_D=0.08)', 'alpha(Z)');
    hold off
end

disp('');

