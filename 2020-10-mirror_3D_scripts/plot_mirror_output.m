function plot_mirror_output(nodes, outputs, pde)

      x = nodes{1};
      x_E = 0.5*3.3443*10^-27*x.^2/(1.602*10^-19);
      num_steps = length(outputs.time_array);
      f1d_analytic = outputs.f_realspace_analytic_nD_t{1,num_steps};
      f1d_ic = outputs.f_realspace_analytic_nD_t{1,1};
for j = 1:num_steps
    f1d = outputs.f_realspace_nD_t{1,j};
     for i = 1:length(x)
         f1d_new(i) = sqrt(x_E(i))*f1d(i);
         f1d_analytic_new(i) = sqrt(x_E(i))*f1d_analytic(i);
     end
     loglog(x_E,f1d_new,'-','LineWidth', 3);
     hold on;
end
     %plot(x,f1d_ic,'--');
     loglog(x_E,f1d_analytic_new, '-o');
     hold off
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
