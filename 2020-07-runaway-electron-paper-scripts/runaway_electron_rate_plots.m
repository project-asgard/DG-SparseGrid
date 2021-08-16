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
T_e_eV = 1e3;
T_e_K = T_e_eV * 11605;
coulomb_ln = 20;

% calculate Dreicer E field in cgs

cgs.n_e = 5e13;
cgs.ln_hat_ee = coulomb_ln;%23.5-log(cgs.n_e^(1/2)*cgs.T_e_eV^(-5/4))-sqrt(10^-5+((log(cgs.T_e_eV)-2)^2)/16);
cgs.th_e = cgs.m_e * cgs.c^2 / (T_e_eV);
% cgs.v_th = 4.19e7*T_e_eV^(1/2);
cgs.v_th = sqrt(cgs.kB*T_e_K/cgs.m_e);
cgs.E_D = 4*pi * cgs.e^3 * cgs.n_e * cgs.ln_hat_ee / (cgs.m_e * cgs.v_th^2); % from Connor Hastie 1975
cgs.E_R = 4*pi*cgs.n_e*cgs.e^3*cgs.ln_hat_ee/(cgs.m_e*cgs.c^2);
cgs.v_c = @(E_cgs) sqrt(4*pi*cgs.e^3 * cgs.n_e * cgs.ln_hat_ee / (cgs.m_e * E_cgs));
disp(['E_D V/m (from cgs): ',num2str(cgs.E_D*300*100)]);
disp(['E_R V/m (from cgs): ',num2str(cgs.E_R*300*100)]);
disp(['E_R V/m (from cgs 2): ',num2str(cgs.E_D*300*100*cgs.kB*T_e_K/(cgs.m_e*cgs.c^2))]);
disp(['E_D / E_R : ', num2str(cgs.E_D/cgs.E_R)]);

% calculate Dreicer E field in SI

SI.n_e = cgs.n_e*1e6;
SI.v_th = sqrt(T_e_eV*SI.e/SI.m_e);
SI.ln_hat_ee = coulomb_ln;
SI.E_D = 1/(4*pi*SI.eps0^2) * SI.n_e * SI.e^3 * SI.ln_hat_ee / (SI.m_e * SI.v_th^2);
disp(['E_D V/m (from SI): ', num2str(SI.E_D)]);

% calculate the RE production rate

cgs.nu_ee2 = 4*pi * cgs.n_e * cgs.e^4 * cgs.ln_hat_ee / (cgs.th_e^(3/2) * cgs.m_e^3 * cgs.c^3);
cgs.nu_ee = 4*pi * cgs.n_e * cgs.e^4 * cgs.ln_hat_ee / (cgs.m_e^2 * cgs.v_th^3);
disp(['nu_ee 1/s (from cgs): ',num2str(cgs.nu_ee)]);
disp(['nu_ee 1/s (from cgs 2): ',num2str(cgs.nu_ee2)]);

% normalized parameters for use in FokkerPlanck2_complete_div

norm.delta = SI.v_th / SI.c;
norm.T_tilde = T_e_eV * SI.e;
norm.v_th_tilde = sqrt(2*norm.T_tilde/SI.m_e);
norm.p = SI.m_e * norm.v_th_tilde;
norm.n_tilde = SI.n_e;
norm.nu_ee_tilde = SI.e^4 * norm.n_tilde * SI.ln_hat_ee / (4*pi*SI.eps0^2*SI.m_e^2*norm.v_th_tilde^3);
norm.tau = norm.nu_ee_tilde;
norm.E_D_tilde = SI.e^3 * norm.n_tilde * SI.ln_hat_ee / (4*pi*SI.eps0^2 * norm.T_tilde);
norm.E = norm.E_D_tilde/2;

disp(['delta: ', num2str(norm.delta)]);
disp(['tau: ', num2str(norm.tau)]); % should this be self-consistent with things in fokkerplanck_paramaters.m?
disp(['E_D: ', num2str(norm.E_D_tilde)]);
disp(['v_th_norm: ', num2str(norm.v_th_tilde * SI.m_e / norm.p)]);

Cn = 1/SI.n_e;

% Franz
T1 = @(E,Z_b) Cn .* cgs.nu_ee .* (E ./ cgs.E_D).^(3/16*Z_b+1);
T2 = @(E,Z_b) exp(-cgs.E_D./(4.*E)-sqrt((1+Z_b).*cgs.E_D./E));

% Connor & Hastie
T1 = @(E,Z_b) Cn .* cgs.n_e * cgs.nu_ee * (E/cgs.E_D)^(-3/8);
T2 = @(E,Z_b) exp(-cgs.E_D/(4*E)-sqrt(2*cgs.E_D/E));

Snr = @(E,Z_b) T1(E,Z_b) .* T2(E,Z_b);

% Scan over E/E_D (0 to 0.6) for Zb=1 - compare with Fig 9 (pg 68) of Franz

ratio_max = 0.3;
N = 100;
ratio = linspace(0.,ratio_max,N);
Z_b = 1;
for i=1:N
   res(i) = Snr(ratio(i)*cgs.E_D,Z_b);
end

figure
semilogy(ratio,res);
ylim([1e-12 1]);


N2 = 15;
ratio2 = linspace(0.,ratio_max,N2);
for i=1:N2
   args.E = ratio2(i)*2; % E into asgard is 2*E/E_D per the normalization
   E_cgs = args.E / 300 / 100;
   v_c_cgs = cgs.v_c(E_cgs+1e-6);
   v_c_SI = v_c_cgs/100;
   v_c_norm = v_c_SI * SI.m_e / norm.p;
   disp(['v_c (cm/s): ',num2str(v_c_cgs)]);
   disp(['v_c (m/s): ',num2str(v_c_SI)]);
   disp(['v_c (norm): ', num2str(v_c_norm)]);
   args.p_max = v_c_norm;

   disp(i);
   [~,~,~,~,~,outputs(i)] = asgard(@fokkerplanck2_complete_div,'timestep_method','BE','num_steps',10,'dt',1,'deg',6,'lev',5,'case',5,'cmd_args',args,'quiet',true,'calculate_mass',true,'grid_type','SG');
   for t=1:numel(outputs(i).time_array)
       pgrid = outputs(i).nodes_t{t}{1};
       zgrid = outputs(i).nodes_t{t}{2};
       f = outputs(i).f_realspace_nD_t{t};
       p2d = repmat(pgrid',1,numel(f(1,:)))';
       p_cutoff = floor(numel(pgrid)/2);
       
       coord = get_realspace_coords(outputs(i).pde,outputs(i).nodes_t{t});
       M1(i,t) = calculate_mass(outputs(i).pde,outputs(i).opts,coord,f(:));
       f(:,1:p_cutoff)=0;
       M2(i,t) = calculate_mass(outputs(i).pde,outputs(i).opts,coord,f(:));
       %M1(i,t) = trapz(zgrid,trapz(pgrid,p2d.^2 .* f,2));
       %M2(i,t) = trapz(zgrid,trapz(pgrid(p_cutoff:end),p2d(:,p_cutoff:end).^2 .* f(:,p_cutoff:end),2));
   end 
%    figure
%    plot(outputs(i).time_array,M1(i,:));
%    hold on
%    plot(outputs(i).time_array,M2(i,:));
%    hold off
%    legend('M1(t)','M2(t)');
%    
%    figure
%    semilogy(outputs(i).time_array,abs(1-M1(i,:)))
%    hold on
%    semilogy(outputs(i).time_array,M2(i,:))
%    hold off
%    legend('1-M1','M2');
%    title('mass error');
   
end

figure
semilogy(ratio,res);
ylim([1e-12 1]);
hold on

norm_fac = res(N)/M2(N2,end);
semilogy(ratio2,M2(:,end)*norm_fac);

M3 = 1-M1;
norm_fac2 = res(N)/M3(N2,end);
semilogy(ratio2,M3(:,end)*norm_fac2);

legend('analytic','M2','1-M1');
hold off


% Scan over Zb (1 to 10) for E=0.08*E_D - compare with Fig 11 (pg 69) of
% Franz

disp('');

