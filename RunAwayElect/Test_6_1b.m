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

% load PDE
PDE_Test_6_1b

FullModel2D_C;

Mat_All = Mat_C;

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


Mat = Inv*Mat_All;
InvMat = inv(speye(DoFs,DoFs)-dt*Mat);


[quad_x,quad_w]=lgwt(num_plot,-1,1);
quad_w = 2^(-Lev)/2*quad_w;
ww = repmat(quad_w,2^Lev,1);
ww = kron(ww,ww)*(LEnd-LInt)*(pEnd-pInt);

val0 = reshape(MM*(F0),nz,nz);
figure
for T = 1 : MaxT

    
    Fn = InvMat*F0;
    F0 = Fn;
    
    time = time + dt;
    
    val_plot = reshape(MM*(F0),nz,nz);
    title(['time = ', num2str(time)])

    surf(x_2D_plot,y_2D_plot,val_plot,val_plot)
    view(0,90);
    colormap jet
    shading interp
    pause(0.1)
    
    val = MM*F0;
    TolTime(T) = time;
    TolMass(T) = sum(ww.*val(:).*y_2D_plot(:).^2);

end

figure;
plot(TolTime,TolMass)
