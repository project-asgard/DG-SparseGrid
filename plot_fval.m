function plot_fval(pde,nodes,f_nD,f_nD_analytic)

num_dims = numel(pde.dimensions);

overPlotAnalytic = 0;
if nargin >= 4
    overPlotAnalytic = 1;
end

if num_dims==1
    
    %%
    % Plot solution
    
    x = nodes{1};
    plot(x,f_nD,'-o');
    
    %%
    % Overplot analytic solution
    
    if pde.checkAnalytic
        hold on;
        coord = {x};
        if overPlotAnalytic
            plot(x,f_nD_analytic,'-');
        end
        hold off;
    end
        
end

if num_dims==2
   
    x = nodes{1};
    y = nodes{2};
    
    nx = numel(x);
    ny = numel(y);
    
    %%
    % Plot a 1D line through the solution
    
    sy = max(1,floor(ny/2));
    if ny > 2
        sy = sy+2; % just to get off the exact middle
    end
    
    f_slice = f_nD(sy,:);
    x = nodes{1};
    y = nodes{2};
    ax1 = subplot(2,2,1);
    plot(x,f_slice,'-o');
    title('1D slice (vertical)');
    
    %%
    % Overplot analytic solution
    
    if pde.checkAnalytic
        f_slice_analytic = f_nD_analytic(sy,:);
        hold on;
        plot(x,f_slice_analytic,'-');
        hold off;
    end
    
    sx = max(1,floor(nx/2));
    if nx > 2
        sx = sx+2; % just to get off the exact middle
    end
    
    f_slice = f_nD(:,sx);
    x = nodes{1};
    y = nodes{2};
    ax1 = subplot(2,2,2);
    plot(y,f_slice,'-o');
    title('1D slice (horizontal)');
    
    if pde.checkAnalytic
        f_slice_analytic = f_nD_analytic(:,sx);
        hold on;
        plot(y,f_slice_analytic,'-');
        hold off;
    end
    
    %%
    % Plot 2D
    
    ax1 = subplot(2,2,3);
    f_nD_with_noise = f_nD;
    f_nD_with_noise(1,1) = f_nD_with_noise(1,1)*1.0001;
    contourf(x,y,f_nD_with_noise,'LineColor','none');
    title('numeric 2D solution');
    
    coordinates = get_realspace_coords(pde,nodes);
    if nargin >= 6
        hold on
        scatter(coordinates(:,1),coordinates(:,2),60,'+','MarkerEdgeColor','white')
        hold off
    end
    
    if pde.checkAnalytic && norm(f_nD_analytic-f_nD_analytic(1,1))>0
        ax2 = subplot(2,2,4);
        contourf(x,y,f_nD_analytic);
        title('analytic 2D solution');
    end
    
%     figure(9)
%     clf
%     [xx,yy] = meshgrid(x,y);
%     p_par = xx.*yy;
%     p_pen = (abs(1-yy.^2)).^(1/2).*xx;
%     f_nD = f_nD_with_noise;
%     contour(p_par,p_pen,f_nD,10,'LineWidth',2)
end

if num_dims==3
    
    figure(1000);
    
    %%
    % Plot a 1D line through the solution
    
    sz = numel(f_nD(:,1,1))/2;
    sy = numel(f_nD(1,:,1))/2;
    sx = numel(f_nD(1,1,:))/2;
    
    f_slice = f_nD(:,sy,sx);
    x = nodes{1};
    y = nodes{2};
    z = nodes{3};
    ax1 = subplot(3,3,1);
    plot(z,f_slice,'-o');
    title('1D slice through 3D');
    
    %%
    % Overplot analytic solution
    
    if pde.checkAnalytic
        f_slice_analytic = f_nD_analytic(:,sy,sx);
        hold on;
        plot(z,f1d_analytic,'-');
        hold off;
    end
    
    %%
    % Plot a 2D xy plane
    
    ax1 = subplot(3,3,4);
    contourf(z,y,f_nD(:,:,sx)');
    title('2D slice through 3D numeric');
    
    if pde.checkAnalytic
        ax2 = subplot(3,3,7);
        contourf(z,y,f_nD_analytic(:,:,sx)');
        title('2D slice through 3D analytic');
    end
    
    %%
    % Plot a 2D xz plane
    
    ax3 = subplot(3,3,5);
    contourf(z,x,squeeze(f_nD(:,sy,:))');
    title('2D slice through 3D numeric');
    
    if pde.checkAnalytic
        ax3 = subplot(3,3,8);
        contourf(z,x,squeeze(f_nD_analytic(:,sy,:))');
        title('2D slice through 3D analytic');
    end
    
    %%
    % Plot a 2D yz plane
    
    ax3 = subplot(3,3,6);
    contourf(y,x,squeeze(f_nD(sz,:,:))');
    title('2D slice through 3D numeric');
    
    if pde.checkAnalytic
        ax3 = subplot(3,3,9);
        contourf(y,x,squeeze(f_nD_analytic(sz,:,:))');
        title('2D slice through 3D analytic');
    end
    
    
end

pause (0.01)

end
