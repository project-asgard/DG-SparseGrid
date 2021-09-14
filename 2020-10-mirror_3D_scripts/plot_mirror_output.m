function plot_mirror_output(nodes, outputs, pde, opts)

      params = mirror_parameters();
      x = nodes{1};
      y = nodes{2};
      nx = numel(x);
      ny = numel(y);
      sx = max(1, floor(nx/2));
      sy = max(1, floor(2*ny/8));
      energy_func_v = @(x) 0.5*pde.params.a.m.*x.^2/(1.602*10^-19);
      energy_func_z = @(z) z.*0 + 2*pi;
      current_func_v = @(x) 2*pi*pde.params.a.Z.*pde.params.e.*x;
      current_func_z = @(x) cos(x);
      coord = get_realspace_coords(pde,nodes);
      mass_func = @(x) x.*0 + sqrt(2*pi);
      zero_func = @(x) x.*0;
      num_dims = numel(pde.dimensions);
      moment_func_nD{1} = energy_func_v;
      moment_func_nD{2} = energy_func_z;
      mass_func_nD{1} = mass_func;
      mass_func_nD{2} = mass_func;
      current_func_nD{1} = current_func_v;
      current_func_nD{2} = current_func_z;
      x_E = 0.5.*params.a.m*x.^2/params.e;
      num_steps = length(outputs.time_array);
      spitzer_conduct = (3/(4*sqrt(2*pi)))*(4*pi*pde.params.eps0)^2*(1.38e-23*pde.params.a.T_eV*11604)^(3/2)...
          /(pde.params.b2.Z*pde.params.e^2*pde.params.a.m^(1/2)*pde.params.ln_delt);
      spitzer_conduct2 = 4*pi*pde.params.eps0^2*(pde.params.a.m*pde.params.a.vth^2)^(3/2) ... 
          /(pde.params.a.m^(1/2)*pde.params.e^2*pde.params.ln_delt*pde.params.b2.Z);
      f1d_analytic = outputs.f_realspace_analytic_nD_t{1,num_steps};
      f1d_ic = outputs.f_realspace_analytic_nD_t{1,1};
      %getting energy density
      lev_vec = [opts.lev, opts.lev];
     % data = struct2cell(outputs);
     % f_data = cell2mat(data{6,1});
     % f_data = reshape(f_data, [num_steps size(f1d_ic)]);
%      [x_grid, y_grid] = meshgrid(outputs.time_array,x_E);
      outputs.time_array(1) = 0;
for i = 1:numel(outputs.time_array)
      f_data(i,:) = outputs.f_realspace_nD_t{i}(ny/4-5,:).*x.^2;
end
      %contourf(outputs.time_array,x_E,real(f_data)',50);
%      f_data = permute(f_data,[1 3 2]);
      
%      contour(outputs.time_array',x',squeeze(f_data(:,:,sy))');
for j = 1:num_steps
    %plot(x,f1d);
    hold on
    fval_realspace = reshape(outputs.f_realspace_nD_t{j}, [numel(f1d_ic) 1]);
    mass_vals(j) = moment_integral(lev_vec, opts.deg, coord, fval_realspace, mass_func_nD, pde.dimensions);
    mass = moment_integral(lev_vec, opts.deg, coord, fval_realspace, mass_func_nD, pde.dimensions);
    conduct_vals(j) = moment_integral(lev_vec,opts.deg,coord,fval_realspace, current_func_nD,pde.dimensions)/(pde.params.E);
    conduct_vals(j) = conduct_vals(j)/spitzer_conduct;
    energy_vals(j) = moment_integral(lev_vec, opts.deg, coord, fval_realspace, moment_func_nD, pde.dimensions)./mass;
  %        vel_vals(i) = sum(outputs.f_realspace_nD_t{1,i}(:).*coord{:,:}*1e-4);
%    x_hint_t = @(t) -x_hint(t,pde.params.a.vel_vals(j),pde.params.a,pde.params.b);
     %e_hint = integral(x_hint_t, 0, outputs.time_array(j)); 
    %hint_func(j) = energy_vals(1)*exp(-e_hint);
     for i = 1:length(x)
          f1d_new(i) = fval_realspace(i);
          f1d_analytic_new(i) = f1d_analytic(i);
      end
     loglog(x_E,f_data,'-','LineWidth', 3);
     hold on
end
      %formulating the Hinton solution
      timespan = [0 1e-3]; %time span in seconds
      x0 = 3e3; %initial energy in eV
      y = @(x) sqrt(1.602e-19.*x./(0.5.*3.3443e-27)); %velocity function in terms of energy
      [t,E] = ode45(@(t,E) (-params.nu_E(y(E),params.a,params.b)-params.nu_E(y(E),params.a,params.b2)).*E,timespan,x0);
     %plot(x,f1d_ic,'--');
     %plotting Hinton solution alongside numerical values
     %x = x./x0;
     %hold on
     %plot(t,E,'r','LineWidth', 2);
     %plot(outputs.time_array,energy_vals,'-o','LineWidth',1);
     %figure
     %plot(outputs.time_array,mass_vals,'-o','LineWidth',2);
     %hold off
     %f_slice = outputs.f_realspace_nD_t{1,num_steps};
     %v_norm = x./pde.params.b.vth;
     %plot(v_norm,f_slice);
     %figure
     %plot(outputs.time_array,conduct_vals,'-','LineWidth',2);
    % hold on
     %loglog(x_E,f1d_analytic_new, '-o');
     %semilogy(outputs.time_array,hint_func,'k');
    % hold off
%      hold off;
%     y = nodes{2};
%     z = nodes{3};
%      nx = numel(x);
%      ny = numel(y);
%      ax1 = subplot(2,2,1);
%      contourf(y,x,squeeze(f_nD{1,1}(:,:))');
%      ax2 = subplot(2,2,2);
%      contourf(y,x,squeeze(f_nD{1,num_steps}(:,:))');
     %title('Pitch vs Velocity');
     
%      sy = numel(f_nD{1,1}(:,1))/4;
%      sx = 1;
%      f1d = f_nD{1,num_steps}(:,sx);
%      f1d_ic = f_nD{1,1}(:,sx);
%      ax3 = subplot(2,2,3);
%      plot(y,f1d,'-o');
%      hold on
%      plot(y,f1d_ic, '-');
%      hold off
%      f1d = f_nD{1,num_steps}(sy,:);
%      f1d_ic = f_nD{1,1}(sy,:);
%      ax4 = subplot(2,2,4);
%      plot(x,f1d,'-o');
%      hold on
%      plot(x,f1d_ic, '-');
     
%     nz = numel(z);     
%     num_dims = 3;
%     subset_dimensions = 2;
%     deg = 6;
%     space_func = @(x) x.*0 + 1;
%     energy_func = @(x) x.^2;
%     perp_func = @(x) sin(x).^2;
%     par_func = @(x) cos(x).^2;
%     coord = get_realspace_coords(pde,nodes);
%     moment_func_nD = {energy_func,perp_func,space_func};
%    % v_perp_temp = moment_integral(pde.get_lev_vec,deg,coord,f_nD,moment_func_nD,pde.dimensions,subset_dimensions);
%     moment_func_nD = {energy_func,par_func,space_func};
%    % v_par_temp = moment_integral(pde.get_lev_vec,deg,coord,f_nD,moment_func_nD,pde.dimensions,subset_dimensions);
%     %%
%     % Plot a 1D line through the solution
%     for i = 1:nx
%         for j = 1:ny
%             v_par(i,j) = x(i).*cos(y(j));
%             v_perp(i,j) = x(i).*sin(y(j));
%         end
%     end
%     
%     sz = numel(f_nD{1,1}(:,1,1))/2 + 1;

    
    %contourf(v_par,v_perp,f_nD{1,1}(:,:,sx));
%     hold on
%     plot(z,v_perp_temp, 'r');
%     plot(z,v_par_temp, 'b');
%     hold off
    %plotting 
%     f1d = f1d(1,:);
%     x = nodes{1};
%     y = nodes{2};
%     z = nodes{3};
%     x = 9.109*10^-31*x.^2/(1.602*10^-19);
    %xlim([0.01,1e4]);
    %ylim([0.01,1e4]);
%     title('Mirror 2D Relaxation');
%     
%     %%
%     
%      %%
%     % Overplot analytic solution
%     
% %     if pde.checkAnalytic
% %         f_1d_analytic = f_nD_analytic(:,sy,sx);
% %         f_1d_analytic = f_1d_analytic(1,:);
%           hold on;
%           f_1d_analytic = f_nD_analytic{1,2};
%           plot(x,f_1d_analytic,'-');
%           hold off;
% %     end
%     
%     sx = max(1,floor(nx/2));
%     if nx > 2
%         sx = sx+2; % just to get off the exact middle
%     end
%     
%     f_slice = f_nD(:,sx);
%     x = nodes{1};
%     y = nodes{2};
%     ax1 = subplot(2,2,2);
%     plot(y,f_slice,'-o');
%     title('1D slice (horizontal)');
%     
%     
%     %%
%     % Plot 2D
%     
%     ax1 = subplot(2,3,1);
%     contourf(z,y,squeeze(f_nD{1,num_steps}(:,:,sx))');
%     title('Space vs. Pitch');
%     ax2 = subplot(2,3,2);
%     contourf(z,x,squeeze(f_nD{1,num_steps}(:,sy,:))');
%     title('Space vs. Velocity');
%     %sx = 46;
%     ax3 = subplot(2,3,3);
% 
%     
%     n_total = zeros(nz, num_steps);
%     for i = 1:num_steps
%          n_total(:,i) = squeeze(sum(f_nD{1,i}, [2 3]));
%     end
%     ax4 = subplot(2,3,5);
%     contourf(z,time_array,squeeze(n_total(:,:))');
%     title('Number Density vs. Time');
    
end
