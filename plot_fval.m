function plot_fval(pde,nodes,fval_realspace,fval_realspace_analytic)

nDims = numel(pde.dimensions);

overPlotAnalytic = 0;
if nargin == 4
    overPlotAnalytic = 1;
end

if nDims==1
    
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

if nDims==2
    
%     figure(1000)
    
    dimensions = pde.dimensions;
    
    deg1=pde.deg;
    lev1=dimensions{1}.lev;
    deg2=pde.deg;
    lev2=dimensions{2}.lev;
    
    dof1=deg1*2^lev1;
    dof2=deg2*2^lev2;
    
    dofD = dof1*dof2;
    assert(dofD==numel(fval_realspace));
    
    %%
    % Reshape dimension ordering is likely a result of kron product dimension
    % ordering. 
    
    f2d = reshape(fval_realspace,dof2,dof1);
    f2d_analytic = reshape(fval_realspace_analytic,dof2,dof1);
    
    x = nodes{1};
    y = nodes{2};
    
    nx = numel(x);
    ny = numel(y);
    
    %%
    % Plot a 1D line through the solution
    
    sy = ny/2;
    
    f1d = f2d(sy,:);
    x = nodes{1};
    y = nodes{2};
    ax1 = subplot(2,2,1);
    semilogy(x,f1d,'-o');
    title('1D slice (vertical)');
    
    %%
    % Overplot analytic solution
    
    if pde.checkAnalytic
        f1d_analytic = f2d_analytic(sy,:);
        hold on;
        plot(x,f1d_analytic,'-');
        hold off;
    end
    
    sx = nx/2;
    
    f1d = f2d(:,sx);
    x = nodes{1};
    y = nodes{2};
    ax1 = subplot(2,2,2);
    plot(y,f1d,'-o');
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
    contourf(x,y,f2d_with_noise);
    title('numeric 2D solution');
    
    if pde.checkAnalytic && norm(f2d_analytic-f2d_analytic(1,1))>0
        ax2 = subplot(2,2,4);
        contourf(x,y,f2d_analytic);
        title('analytic 2D solution');
    end
    
end

if nDims==3
    
    figure(1000);
    
    dimensions = pde.dimensions;
    
    deg1=pde.deg;
    lev1=dimensions{1}.lev;
    deg2=pde.deg;
    lev2=dimensions{2}.lev;
    deg3=pde.deg;
    lev3=dimensions{3}.lev;
    
    dof1=deg1*2^lev1;
    dof2=deg2*2^lev2;
    dof3=deg3*2^lev3;
    
    dofD = dof1*dof2*dof3;
    assert(dofD==numel(fval_realspace));
    
    f3d = reshape(fval_realspace,dof3,dof2,dof1);
    f3d_analytic = reshape(fval_realspace_analytic,dof3,dof2,dof1);
    
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