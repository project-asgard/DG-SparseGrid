function plot_fval(pde,nodes,f_nD,f_nD_analytic,element_coordinates,window_title)

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
    
    if ~isempty(pde.solutions)
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
%     sy=ny;
    
    f_slice = f_nD(sy,:);
    x = nodes{1};
    y = nodes{2};
    ax1 = subplot(2,2,1);
    plot(x,f_slice,'-o');
    title('1D slice (horizontal)');
    
    %%
    % Overplot analytic solution
    
    if ~isempty(pde.solutions)
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
    title('1D slice (vertical)');
    
    if ~isempty(pde.solutions)
        f_slice_analytic = f_nD_analytic(:,sx);
        hold on;
        plot(y,f_slice_analytic,'-');
        hold off;
    end
    
    %%
    % Plot 2D
    
    ax1 = subplot(2,2,3);
    contourf(x,y,f_nD,'LineColor','none');
    hold off

    title('numeric 2D solution');
    
    if nargin >= 5
        hold on
        scatter(element_coordinates(:,1),element_coordinates(:,2),60,'+','MarkerEdgeColor','black')
        hold off
    end
    
    if ~isempty(pde.solutions) && norm(f_nD_analytic-f_nD_analytic(1,1))>0
        ax2 = subplot(2,2,4);
        contourf(x,y,f_nD_analytic);
        title('analytic 2D solution');
    end

    do_RE_paper_plots = false;
    if do_RE_paper_plots
        figure(11)
        levs = [-9,-8,-7,-6,-5,-4,-3,-2,-1,0];
        f_nD(f_nD<1e-12)=1e-12;
        [M,c]=contourf(x,y,log10(f_nD),levs,'LineColor','none');
        xlabel('p');
        ylabel('\zeta');
        clabel(M,c,levs,'Color','w');
        colormap(flipud(pink));
        set(gca,'FontSize',16)
        g=100;
        yy=yline(0.984,'Color','b','LineWidth',1);
        xx=xline(5.01,'Color','b','LineWidth',1);
        hold on
        scatter(element_coordinates(:,1),element_coordinates(:,2),60,'+','MarkerEdgeColor','r')
        hold off
    end
    
    plot_fval_in_cyl = false;
    if plot_fval_in_cyl
        p = x;
        z = y;
        f = f_nD;
        pper = linspace(0,max(p),31);
        ppar = linspace(-max(p),+max(p),51);
        [ppar2d,pper2d] = meshgrid(ppar,pper);
        p2dA = sqrt(ppar2d.^2+pper2d.^2);
        z2dA = cos(atan2(pper2d,ppar2d));
        f2d = interp2(p,z,f,p2dA,z2dA,'linear',1e-12);
        figure(87)
        f2d(f2d<1e-12)=1e-12;
        contourf(ppar,pper,log10(f2d),levs,'LineColor','none')
        colormap(flipud(pink));
        sg_th = acos(element_coordinates(:,2));
        sg_pper = element_coordinates(:,1).*sin(sg_th);
        sg_ppar = element_coordinates(:,1).*cos(sg_th);
        hold on
        scatter(sg_ppar,sg_pper,60,'+','MarkerEdgeColor','k','LineWidth',2);
        nL = 100;
        pL1 = linspace(0,1,nL)*(max(x)-min(x))+min(x);
        zL1 = zeros(1,nL)+0.984;
        pperL1 = pL1.*sin(acos(zL1));
        pparL1 = pL1.*cos(acos(zL1));
        plot(pparL1,pperL1,'-r','LineWidth',1);
        pL2 = zeros(1,nL)+5.01;
        zL2 = linspace(0,1,nL)*(max(y)-min(y))+min(y);
        pperL2 = pL2.*sin(acos(zL2));
        pparL2 = pL2.*cos(acos(zL2));
        plot(pparL2,pperL2,'-b','LineWidth',1);
        
        pL2 = zeros(1,nL)+max(x);
        zL2 = linspace(0,1,nL)*(max(y)-min(y))+min(y);
        pperL2 = pL2.*sin(acos(zL2));
        pparL2 = pL2.*cos(acos(zL2));
        plot(pparL2,pperL2,'-k','LineWidth',0.5);
        pL2 = zeros(1,nL)+min(x);
        zL2 = linspace(0,1,nL)*(max(y)-min(y))+min(y);
        pperL2 = pL2.*sin(acos(zL2));
        pparL2 = pL2.*cos(acos(zL2));
        plot(pparL2,pperL2,'-k','LineWidth',0.5);
        hold off

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
    
%    figure(1000);
    
    %%
    % Plot a 1D line through the solution
    
    sz = numel(f_nD(1,1,:))/2;
    sy = numel(f_nD(1,:,1))/2;
    sx = numel(f_nD(:,1,1))/2;
%     
%     %plotting in x-direction
%     f1d = f_nD(sz,sy,:);
%     f1d = f1d(1,:);
    x = nodes{1};
    y = nodes{2};
    z = nodes{3};
%     ax1 = subplot(3,3,1);
%     plot(x,f1d,'-o');
%     title('1D slice through velocity dimension');
%     
%     %plotting y-direction
%     f1d = f_nD(sz,:,sx);
%     x = nodes{1};
%     y = nodes{2};
%     z = nodes{3};
%     ax2 = subplot(3,3,2);
%     plot(y,f1d,'-o');
%     title('1D slice through pitch dimension');
%     
%     %plotting z-direction
%     f1d = f_nD(:,sy,sx);
%     x = nodes{1};
%     y = nodes{2};
%     z = nodes{3};
%     ax3 = subplot(3,3,3);
%     plot(z,f1d,'-o');
%     title('1D slice through spatial dimension');
%     
    %%
    % Overplot analytic solution
    
    if ~isempty(pde.solutions)
        f_slice_analytic = f_nD_analytic(:,sy,sx);
        hold on;
        plot(z,f_slice_analytic,'-');
        hold off;
    end
    
    %%
    % Plot a 2D xy plane
    
    ax1 = subplot(2,3,1);
    contourf(z,y,f_nD(:,:,sx)');
    title('2D slice through 3D numeric');
    
    if ~isempty(pde.solutions)
        ax2 = subplot(2,3,4);
        contourf(z,y,f_nD_analytic(:,:,sx)');
        title('2D slice through 3D analytic');
    end
    
    %%
    % Plot a 2D xz plane
    
    ax3 = subplot(2,3,2);
    contourf(z,x,squeeze(f_nD(:,sy,:))');
    title('2D slice through 3D numeric');
    
    if ~isempty(pde.solutions)
        ax3 = subplot(2,3,5);
        contourf(z,x,squeeze(f_nD_analytic(:,sy,:))');
        title('2D slice through 3D analytic');
    end
    
    %%
    % Plot a 2D yz plane
    
    ax3 = subplot(2,3,3);
    contourf(y,x,squeeze(f_nD(sz,:,:))');
    title('2D slice through 3D numeric');
    
    if ~isempty(pde.solutions)
        ax3 = subplot(2,3,6);
        contourf(y,x,squeeze(f_nD_analytic(sz,:,:))');
        title('2D slice through 3D analytic');
    end
    
    
end

pause (0.01)

end
