function plot_fval(pde,nodes,fval_realspace,fval_realspace_analytic,Meval,coordinates)

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
    
%     dof1=deg1*2^lev1;
%     dof2=deg2*2^lev2;
    dof1=numel(Meval{1}(:,1));
    dof2=numel(Meval{2}(:,1));
    
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
    
    stat = mkdir('output_adapt');
    fName = ['output_adapt/f2d_adapt-' sprintf('deg_%d',deg1) '.mat'];
    save(fName, 'x', 'y', 'f2d');
    %%
    % Plot a 1D line through the solution
    
    sy = max(1,ceil(ny-1));
   % if ny > 2
   % sy = sy+2; % just to get off the exact middle
   % end
    
    f1d = f2d(sy,:);
    x = nodes{1};
    y = nodes{2};
    ax1 = subplot(2,2,1);
    %plot(x,f1d,'-o');
    %load('f1d_FG_upper_1,5_deg4', 'x_FG_deg4', 'f1d_FG_deg4');
    %load('f1d_FG_upper_1,5_deg7', 'x_FG_deg7', 'f1d_FG_deg7');
    fNamex = ['output_adapt/f1d_x-' sprintf('deg_%d',deg1) '.mat'];
    save(fNamex, 'x', 'f1d');
    %load('f1d_FG_upper_1,5_deg5', 'x_FG_deg5', 'f1d_FG_deg5');
    %load('f1d_FG_upper_1,5_ref', 'x_FG_ref', 'f1d_FG_ref');
    %load('f1d_SG_upper_1,5s_unadapted', 'x_SG_unadapt', 'f1d_SG_unadapt');
    %load('f1d_SG_upper_1,5s_ref2', 'x_SG_ref2', 'f1d_SG_ref2');
    %load('f1d_FG_upper_1,5_coarse', 'x_FG_coarse', 'f1d_FG_coarse');
   % load('f1d_SG_upper_1,8s_ref1', 'x_SG_ref1', 'f1d_SG_ref1');
    %load('f1d_SG_upper_1,8s_ref', 'x_SG_ref', 'f1d_SG_ref');
    %load('f1d_FG_upper_1,8s', 'x_FG', 'f1d_FG');
    %load('f1d_SG_upper_1,8s', 'x_SG', 'f1d_SG');
    
   % load('f1d_SG', 'x_SG', 'f1d_SG');
   % load('f1d_SG_mid', 'x_SG_mid','f1d_SG_mid');
   % load('f1d_FG', 'x_FG', 'f1d_FG');
   % load('f1d_FG_mid_neu_minusE2', 'x_FG_mid_neu', 'f1d_FG_mid_neu');
   % semilogy(x_SG, f1d_SG, '--m', 'LineWidth', 2);
    ylim([10^-4, 10^0]);
  %  semilogy(x_FG, f1d_FG, '-m', 'LineWidth', 2);
  %  semilogy(x_FG_ref, f1d_FG_ref, '-k', 'LineWidth', 2);
  %  semilogy(x_FG_deg4, f1d_FG_deg4, '-r', 'LineWidth', 2);
   % hold on;
 %   semilogy(x_FG_coarse, f1d_FG_coarse, '-g', 'LineWidth', 2);
  %  semilogy(x_FG_deg5, f1d_FG_deg5, '-m', 'LineWidth', 2);
  %  semilogy(x_FG_deg6, f1d_FG_deg6, '-k', 'LineWidth', 2);
  %  semilogy(x_FG_deg7, f1d_FG_deg7, '-ob');
   % semilogy(x_FG_deg7, f1d_FG_deg7, '-og');
   % semilogy(x_SG_ref, f1d_SG_ref, 'b--', 'LineWidth', 2);
   % semilogy(x_SG_ref1, f1d_SG_ref1, 'r--', 'LineWidth', 2);
   % semilogy(x_SG_ref2, f1d_SG_ref2, 'k--', 'LineWidth', 2);
   % semilogy(x_SG_unadapt, f1d_SG_unadapt, 'c--', 'LineWidth',2);
   % hold off;
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
    fNamey = ['output_adapt/f1d_y_adapt-' sprintf('deg_%d',deg1) '.mat'];
    save(fNamey, 'y', 'f1d');
    %load('f1d_FG_deg5_Vert', 'y_FG_deg5_vert', 'f1d_FG_deg5_vert');
    %load('f1d_FG_deg4_Vert', 'y_FG_deg4_vert', 'f1d_FG_deg4_vert');
    ax1 = subplot(2,2,2);
   % plot(y_FG_deg4_vert,f1d_FG_deg4_vert,'-r', 'LineWidth', 2);
  %  hold on;
  %  plot(y_FG_deg5_vert,f1d_FG_deg5_vert,'-m', 'LineWidth', 2);
  %  plot(y_FG_deg6_vert,f1d_FG_deg6_vert, '-k', 'LineWidth', 2);
  %  hold off;
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
