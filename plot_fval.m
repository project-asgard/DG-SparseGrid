function plot_fval(pde,nodes,fval_realspace,fval_realspace_analytic,Meval,coordinates)

num_dims = numel(pde.dimensions);

overPlotAnalytic = 0;
if nargin >= 4
    overPlotAnalytic = 1;
end

if num_dims==1
    
    %%
    % Plot solution
    
    f1d = fval_realspace;
    x = nodes{1};
    plot(x,f1d,'-o');
    
    %%
    % Overplot analytic solution
    
    if pde.checkAnalytic
        hold on;
        coord = {x};
        %         f1d_analytic = getAnalyticSolution_D(coord,time,pde);
        if overPlotAnalytic
            plot(x,fval_realspace_analytic,'-');
        end
        hold off;
    end
        
end

if num_dims==2
            
    f2d = singleD_to_multiD(num_dims,fval_realspace,nodes);
    f2d_analytic = singleD_to_multiD(num_dims,fval_realspace_analytic,nodes);
   
    x = nodes{1};
    y = nodes{2};
    
    nx = numel(x);
    ny = numel(y);
    
    %%
    % Plot a 1D line through the solution
    
    sy = max(1,floor(ny-1));
%    if ny > 2
%        sy = sy+2; % just to get off the exact middle
%    end
    
    f1d = f2d(sy,:);
    x = nodes{1};
    y = nodes{2};
    ax1 = subplot(2,2,1);
    plot(x,f1d,'-o');
    semilogy(x,f1d,'LineWidth', 2); %semilog for SG output
    ylim([10^-10, 10^0]);
    title('1D slice (vertical)');
    
    %%
    % Overplot analytic solution
    
    if pde.checkAnalytic
        f1d_analytic = f2d_analytic(sy,:);
        hold on;
        plot(x,f1d_analytic,'-');
        hold off;
    end
    
    sx = max(1,floor(nx/2));
    if nx > 2
        sx = sx+2; % just to get off the exact middle
    end
    
    f1d = f2d(:,sx);
    x = nodes{1};
    y = nodes{2};
    ax1 = subplot(2,2,2);
    plot(y,f1d,'-o');
    semilogy(y,f1d,'LineWidth',2); %semilog for SG output
    title('1D slice (horizontal)');
    
    if pde.checkAnalytic
        f1d_analytic = f2d_analytic(:,sx);
        hold on;
        plot(y,f1d_analytic,'-');
        hold off;
    end
    
    %%
    % Plot 2D
    
    ax1 = subplot(2,2,3);
    f2d_with_noise = f2d;
    f2d_with_noise(1,1) = f2d_with_noise(1,1)*1.0001;
    contourf(x,y,f2d_with_noise,'LineColor','none');
    title('numeric 2D solution');
    xline(x(sx), 'LineWidth', 2);
    yline(y(sy), 'LineWidth', 2);
    if nargin >= 6
        hold on
        scatter(coordinates(:,1),coordinates(:,2),60,'+','MarkerEdgeColor','white')
        hold off
    end
    
    if pde.checkAnalytic && norm(f2d_analytic-f2d_analytic(1,1))>0
        ax2 = subplot(2,2,4);
        contourf(x,y,f2d_analytic);
        title('analytic 2D solution');
    end
    
    figure(9)
    clf
    [xx,yy] = meshgrid(x,y);
    p_par = xx.*yy;
    p_pen = (abs(1-yy.^2)).^(1/2).*xx;
    f2d = f2d_with_noise;
    contour(p_par,p_pen,f2d,10,'LineWidth',2)
end

if num_dims==3
    
    figure(1000);
    
    dimensions = pde.dimensions;
    
    f3d = singleD_to_multiD(num_dims,fval_realspace,nodes);
    f3d_analytic = singleD_to_multiD(num_dims,fval_realspace_analytic,nodes);
    
    %%
    % Plot a 1D line through the solution
    
    sz = numel(f3d(:,1,1))/2;
    sy = numel(f3d(1,:,1))/2;
    sx = numel(f3d(1,1,:))/2;
    
    f1d = f3d(:,sy,sx);
    x = nodes{1};
    y = nodes{2};
    z = nodes{3};
    ax1 = subplot(3,3,1);
    plot(z,f1d,'-o');
    title('1D slice through 3D');
    
    %%
    % Overplot analytic solution
    
    if pde.checkAnalytic
        f1d_analytic = f3d_analytic(:,sy,sx);
        hold on;
        plot(z,f1d_analytic,'-');
        hold off;
    end
    
    %%
    % Plot a 2D xy plane
    
    ax1 = subplot(3,3,4);
    contourf(z,y,f3d(:,:,sx)');
    title('2D slice through 3D numeric');
    
    if pde.checkAnalytic
        ax2 = subplot(3,3,7);
        contourf(z,y,f3d_analytic(:,:,sx)');
        title('2D slice through 3D analytic');
    end
    
    %%
    % Plot a 2D xz plane
    
    ax3 = subplot(3,3,5);
    contourf(z,x,squeeze(f3d(:,sy,:))');
    title('2D slice through 3D numeric');
    
    if pde.checkAnalytic
        ax3 = subplot(3,3,8);
        contourf(z,x,squeeze(f3d_analytic(:,sy,:))');
        title('2D slice through 3D analytic');
    end
    
    %%
    % Plot a 2D yz plane
    
    ax3 = subplot(3,3,6);
    contourf(y,x,squeeze(f3d(sz,:,:))');
    title('2D slice through 3D numeric');
    
    if pde.checkAnalytic
        ax3 = subplot(3,3,9);
        contourf(y,x,squeeze(f3d_analytic(sz,:,:))');
        title('2D slice through 3D analytic');
    end
    
    
end

pause (0.01)

end
