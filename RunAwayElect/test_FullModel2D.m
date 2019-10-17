addpath(genpath(pwd))

Lev = 5;
Deg = 2;
num_plot = 2;

DoFs = (2^Lev*Deg);
close all

LInt = -1;
LEnd = 1;

pInt = 0.1;
pEnd = 10;

Lmax = LEnd-LInt;
dx = Lmax/2^Lev;

CFL = 0.001;
dt = .1;%CFL*(dx)^3;
MaxT = ceil(1e3/dt);


FullModel2;

FullModel2D_C;
FullModel2D_E;
FullModel2D_R;

Mat_All = Mat_C;%Mat_C + Mat_E;% + Mat_R; %

Mat_Mass = kron(speye(DoFs,DoFs),Mat_Mass_p);
Inv = inv(Mat_Mass);


DoFs = (2^Lev*Deg)^2;

F0 = zeros(DoFs,1);

% Initial Condition
time = 0;
F0 = ComputRHS2D(Lev,Deg,[LInt pInt ],[LEnd pEnd ],Exa0,time);

[x_node,Meval] = PlotDGData(Lev,Deg,LInt,LEnd,num_plot);
[p_node,Peval] = PlotDGData(Lev,Deg,pInt,pEnd,num_plot);
[x_2D_plot,y_2D_plot] = meshgrid(x_node,p_node);
MM = kron(Meval,Peval);
nz = size(x_2D_plot,1);

val_plot = reshape(MM*(F0),nz,nz);
valval0 = val_plot;
subplot(3,1,1)
surf(x_2D_plot,y_2D_plot,val_plot,val_plot)
shading interp
view(0,90);

subplot(3,1,2)
semilogy(y_2D_plot(:,4),abs(val_plot(:,4)))
%     axis([0 5 1e-9 1])

subplot(3,1,3)
plot(x_2D_plot(4,:),val_plot(4,:))
pause(0.1)
colormap jet

p_par = x_2D_plot.*y_2D_plot;
p_pen = (1-x_2D_plot.^2).^(1/2).*y_2D_plot;
figure
surf(p_par,p_pen,val_plot,val_plot)
shading interp
view(0,90);

Mat = Inv*Mat_All;
InvMat = inv(speye(DoFs,DoFs)-dt*Mat);


[quad_x,quad_w]=lgwt(num_plot,-1,1);
quad_w = 2^(-Lev)/2*quad_w;
ww = repmat(quad_w,2^Lev,1);
ww = kron(ww,ww)*(LEnd-LInt)*(pEnd-pInt);

val0 = reshape(MM*(F0),nz,nz);
figure
for T = 1 : MaxT
    
    %     val_plot = reshape(MM*(F0),nz,nz);
    %     subplot(2,1,1)
    %     surf(x_2D_plot,y_2D_plot,val_plot,val_plot)
    %     shading interp
    %     colormap jet
    %     view(-90,90);
    %     title(['time = ', num2str(time)])
    %
    %     subplot(2,1,2)
    %     plot(val_plot(:,4))
    %     pause(0.1)
    
    %     F1 = F0 + dt*(Mat)*F0;
    
    % RK3
    %     F1 = F0 + dt*(  Mat*F0 );
    %     F2 = 3/4*F0+1/4*F1+1/4*dt*(Mat*F1);
    %     Fn = 1/3*F0+2/3*F2+2/3*dt*(Mat*F2);
    %     F0 = Fn;
    
    Fn = InvMat*F0;
    F0 = Fn;
    
    time = time + dt;
    
    val_plot = reshape(MM*(F0),nz,nz);
    subplot(3,1,1)
    %     surf(x_2D_plot,y_2D_plot,val_plot-valval0,val_plot-valval0)
    %     shading interp
    %     colormap jet
    %     surf(p_par,p_pen,val_plot,real(log10(val_plot)))
    %         view(0,90);
    %             colormap jet
    % shading interp
    contour(p_par,p_pen,val_plot,20,'LineWidth',2)
    %     shading interp
    %
    %     colormap jet
    %     view(0,90);
    title(['time = ', num2str(time)])
    
    subplot(3,1,2)
    %     plot(y_2D_plot(:,4),val_plot(:,4))
    semilogy(y_2D_plot(:,4),abs(val_plot(:,4)-valval0(:,4)))
    %     axis([0 10 1e-10 1])
    
    subplot(3,1,3)
    %     plot(x_2D_plot(4,:),val_plot(4,:)-valval0(4,:))
    surf(p_par,p_pen,val_plot,val_plot)
    view(0,90);
    colormap jet
    shading interp
    pause(0.1)
    
    val = MM*F0;
    TolTime(T) = time;
        TolMass(T) = sum(ww.*val(:).*y_2D_plot(:).^2);
%     TolMass(T) = sum(ww.*val(:).*p_pen(:));
end

figure;
plot(TolTime,TolMass)
