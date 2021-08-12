% Produce the runaway electron fraction plots

% Calculate the relativistic and non-relativistic runaway electron rates
% from what is in Franz's thesis (pg 71 - which comes from Connor and Hastie, Nucl.
% Fusion, 15:415, 1975).

e = 1.6e-19;
c = 3e8;
m_e = 9.1e-31;
T_e = 1e3*e;
n_e = 5e19;
Z_b = 1;
ln_hat_ee = 15;
th_e = m_e * c^2 / T_e;
v_th = sqrt(T_e/m_e);
E_D = 2*pi * e^3 * n_e * ln_hat_ee / (m_e * v_th^2);
nu_ee = 4*pi * n_e * e^4 * ln_hat_ee / (th_e^(3/2) * m_e^3 * c^3);
Cn = 0.5e-68;

T1 = @(E,Z_b) Cn .* nu_ee .* (E ./ E_D).^(3/16*Z_b+1);
T2 = @(E,Z_b) exp(-E_D./(4.*E)-sqrt((1+Z_b).*E_D./E));

Snr = @(E,Z_b) T1(E,Z_b) .* T2(E,Z_b);

% Scan over E/E_D (0 to 0.6) for Zb=1 - compare with Fig 9 (pg 68) of Franz

N = 100;
ratio = linspace(0.,0.6,N);
Z_b = 1;
for i=1:N
   res(i) = Snr(ratio(i)*E_D,Z_b);
end



N2 = 10;
ratio2 = linspace(0.,0.6,N2);
for i=1:N2
   args.E = ratio2(i)*1;%E_D; 
   disp(i);
   [~,~,~,~,~,outputs(i)] = asgard(@fokkerplanck2_complete_div,'timestep_method','matrix_exponential','num_steps',1,'dt',10.0,'deg',3,'lev',3,'case',5,'cmd_args',args,'quiet',true);
   for t=1:numel(outputs(i).time_array)
       pgrid = outputs(i).nodes_t{t}{1};
       zgrid = outputs(i).nodes_t{t}{2};
       f = outputs(i).f_realspace_nD_t{t};
       p2d = repmat(pgrid',1,numel(f(1,:)))';
       p_cutoff = numel(pgrid)/2;
       M1(i,t) = trapz(zgrid,trapz(pgrid,p2d.^2 .* f,2));
       M2(i,t) = trapz(zgrid,trapz(pgrid(p_cutoff:end),p2d(:,p_cutoff:end).^2 .* f(:,p_cutoff:end),2));
   end   
%    plot(outputs(i).time_array,M1(i,:));
%    hold on
%    plot(outputs(i).time_array,M2(i,:));
%    hold off
end

figure
semilogy(ratio,res*1e60);
ylim([1e-12 1]);
hold on
semilogy(ratio2,M2(:,end)*1e-2);
hold off

% Scan over Zb (1 to 10) for E=0.08*E_D - compare with Fig 11 (pg 69) of
% Franz

disp('');

