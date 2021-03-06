function plot_mirror_output(nodes, f_nD, f_nD_analytic, pde)

    x = nodes{1};
    y = nodes{2};
    z = nodes{3};
    
    nx = numel(x);
    ny = numel(y);
    nz = numel(z);
    
    %%
    % Plot a 1D line through the solution
    
    sz = numel(f_nD(:,1,1))/2;
    sy = numel(f_nD(1,:,1))/2;
    sx = numel(f_nD(1,1,:)) - 1;
    
    %plotting 
    f1d = f_nD(:,sy,sx);
    f1d = f1d(1,:);
    x = nodes{1};
    y = nodes{2};
    z = nodes{3};
    x = 9.109*10^-31*x.^2/(1.602*10^-19);
    plot(z,f1d,'-o');
    %xlim([0.01,1e4]);
    %ylim([0.01,1e4]);
    title('Mirror 2D Relaxation');
    
    %%
    
     %%
    % Overplot analytic solution
    
    if pde.checkAnalytic
        f_1d_analytic = f_nD_analytic(:,sy,sx);
        f_1d_analytic = f_1d_analytic(1,:);
        hold on;
        plot(z,f_1d_analytic,'-');
        hold off;
    end
    
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
%     ax1 = subplot(2,2,3);
%     f_nD_with_noise = f_nD;
%     f_nD_with_noise(1,1) = f_nD_with_noise(1,1)*1.0001;
%     contourf(x,y,f_nD_with_noise,'LineColor','none');
%     title('numeric 2D solution');
%     


end
