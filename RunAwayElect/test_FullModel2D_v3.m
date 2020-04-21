addpath(genpath(pwd))
close all
clc
clear

% We may consider the upwind scheme

Lev = 5;
Deg = 2;
num_Int = Deg+1;
num_plot = Deg+1;

CutOff = 10;

DoFs = (2^Lev*Deg);


LInt = -1;
LEnd = 1;

pInt = 0;
pEnd = 20;

Lmax = LEnd-LInt;
dx = Lmax/2^Lev;

CFL = 0.001;
dt = 1e-3;%1e-3;
% MaxT = ceil(3e2/dt);
MaxT = ceil(2e2/dt);

FullModel2_v2;
FullModel2D_C;
FullModel2D_E;
FullModel2D_R;

Mat_All = Mat_C + Mat_R + Mat_E ;

Mat_Mass = kron(speye(DoFs,DoFs),Mat_Mass_p);
% Inv = inv(Mat_Mass);


DoFs = (2^Lev*Deg)^2;

F0 = zeros(DoFs,1);

% Initial Condition
time = 0;
F0 = ComputRHS2D(Lev,Deg,[LInt pInt ],[LEnd pEnd ],Exa0,time);

[x_node,Meval] = PlotDGData(Lev,Deg,LInt,LEnd,num_plot);
[p_node,Peval] = PlotDGData(Lev,Deg,pInt,pEnd,num_plot);

% quadrature for integration
[x_node2,Meval2] = PlotDGData(Lev,Deg,LInt,LEnd,num_Int);
[p_node2,Peval2] = PlotDGData(Lev,Deg,pInt,pEnd,num_plot);

[x_2D_plot,y_2D_plot] = meshgrid(x_node,p_node);
[x_2D_Int,y_2D_Int] = meshgrid(x_node2,p_node2);
MM = kron(Meval,Peval);

MMInt = kron(Meval2,Peval2);
nz = size(x_2D_plot,1);

val_plot = reshape(MM*(F0),nz,nz);
tmp0 = reshape(val_plot,num_plot*2^Lev,num_plot*2^Lev);
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

% Mat = Inv*Mat_All;
Mat = Mat_Mass\Mat_All;
% InvMat = inv(speye(DoFs,DoFs)-dt*Mat);
InvMat = (speye(DoFs,DoFs)-dt*Mat);



[quad_x,quad_w]=lgwt(num_plot,-1,1);
quad_w = 2^(-Lev)/2*quad_w;
w_x = repmat(quad_w,2^Lev,1)/2^Lev;
ww = repmat(quad_w,2^Lev,1);
ww = kron(ww,ww)*(LEnd-LInt)*(pEnd-pInt);
ww2 = ww;

T = 0;
time = 0;
val = MMInt*F0;
TolTime(T+1) = time;
TolMass(T+1) = sum(ww2.*val(:).*y_2D_Int(:).^2);
iy = y_2D_Int(:)>=CutOff;
    
TolMass2(T+1) = sum(ww2(iy).*val(iy).*y_2D_Int(iy).^2);
    
figure
for T = 1 : MaxT
    
    Fn = InvMat\F0;
    F0 = Fn;
    
    time = time + dt;
    
    val_plot = reshape(MM*(F0),nz,nz);
    
    figure(10)
    surf(x_2D_plot,y_2D_plot,val_plot,val_plot)
    shading interp
    colormap jet
    view(0,90)
    
    figure(2)
    [c,h] = contour(p_par,p_pen,val_plot,10,'LineWidth',2);
    colormap jet
    title(['time = ', num2str(time)])
    
    figure(4)
    %     plot(x_2D_plot(4,:),val_plot(4,:)-valval0(4,:))
    surf(p_par,p_pen,val_plot,val_plot)
    view(0,90);
    colorbar
    colormap jet
    shading interp
%     
    

%     if mod(T,100) == 0
%         save(['time=',num2str(time),'.mat'],'F0')
%     end
    
    val = MMInt*F0;
    TolTime(T+1) = time;
    TolMass(T+1) = sum(ww2.*val(:).*y_2D_Int(:).^2);
    iy = y_2D_Int(:)>=CutOff;
    
    TolMass2(T+1) = sum(ww2(iy).*val(iy).*y_2D_Int(iy).^2);
    
%     figure(3);
%     plot(TolTime,TolMass,'r-','LineWidth',2);hold on;
%     plot(TolTime,TolMass2,'b-','LineWidth',2);hold off;
    
    tmp = reshape(val_plot,num_plot*2^Lev,num_plot*2^Lev);
    
    figure(12)
    semilogy(p_node,abs(tmp*w_x),'b-','LineWidth',2);hold on
    semilogy(p_node,abs(tmp0*w_x),'r-','LineWidth',2);hold off
    axis([0 pEnd 1e-7 3])
    
    pause(0.1)
    
    [T,TolMass2(T+1),TolMass(T+1),time]
    if TolMass(T+1)<0
        break
    end

end

figure;
plot(TolTime,TolMass,'LineWidth',2)
hold on;
plot(TolTime,TolMass2,'LineWidth',2)
