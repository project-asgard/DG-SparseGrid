clear all
close all
% clc

% Test
sigma = 0.1;
f0 = @(x)( exp(-x.^2/sigma^2) );
% f0 = @(x)(x-x+1);
phi = @(x,t)( tanh(atanh(x)-t) );
exactf = @(x,t)(...
    (1-phi(x,t).^2)./(1-x.^2).*f0(phi(x,t)) ...
    );
funcCoef = @(x)(1-x.^2);

format short e
addpath(genpath(pwd))

% start with coarse mesh
Lev = 4;
Deg = 2;
num_plot = Deg;

Lstart = -1;
Lend = 1;
EndTime = 3;

CFL = 0.01;
% dt = CFL*h^((Deg-1)/3)/2;
h = (Lend-Lstart)/2^Lev;
dt = CFL*h^((Deg)/3);
maxT = ceil(EndTime/dt)

load(['Mat_Lev10_Deg',num2str(Deg),'_v2.mat']);

% [Hash,IHash,LeafHash,ILeaf,FineIndex] = HashTable1D(Lev);
[Hash,IHash,FineIndex] = HashTable1D(Lev);
DoFs = Deg*2^Lev;

Mat_tmp = Mat(1:DoFs,1:DoFs);


quad_num = Deg;
[quad_x,quad_w]=lgwt(quad_num,-1,1);
p_val = legendre(quad_x,Deg);
% compute initial condition
for L=0:2^Lev-1
    
    x0 = Lstart+L*h;
    x1 = x0+h;
    xi = quad_x*(x1-x0)/2+(x1+x0)/2;
    
    val = sqrt(h)/2*[p_val'*(quad_w.*exactf(xi,0))];
    
    c = Deg*L+1:Deg*(L+1);
    fval(c,1) = val;

end

for i = 1:size(IHash,2)
    l1 = IHash{i};
    Id(Deg*(i-1)+[1:Deg]) = Deg*(l1(3)-1)+[1:Deg];
end

num_In = 2^5;%2^(10-Lev);
plot_start = 2^(10-Lev)/2;
MMeval = Meval(plot_start:num_In:end,Id);
xx_node = x_node(plot_start:num_In:end);

clear FMWT
FMWT = OperatorTwoScale(Deg,2^Lev);
% Meval2 = Meval2*FMWT(1:DoFs,1:DoFs)';

fval_MW = FMWT*fval;

% construct x_node through Hash table
% % x_node = zeros(size(IHash,2),1);
% % for i = 1:size(IHash,2)
% %         l1 = IHash{i};
% %         key = [l1(1),l1(2)];
% %         lr = Hash.(sprintf('i%g_',key));
% %         
% %         x_node
% %         
% %         Id(Deg*(i-1)+[1:Deg]) = Deg*(lr-1)+[1:Deg];
% % end

% return
figure(1);
% subplot(1,2,1)
% plot(xx_node,MMeval*fval_MW,'r-o',xx_node,exactf(xx_node,0),'b--','LineWidth',2)
% % hold on;
% % plotgrid(IHash,Lstart,Lend,-0.1);
% pause(0.1)

% return
% figure(2);
% plot(x_node2,Meval2*fval_MW,'r-o',x_node2,exactf(x_node2,0),'b--','LineWidth',2)
% hold on;
% plotgrid(IHash,Lstart,Lend,-0.1);
% pause(0.1)
% return
return
LeafSolIndex = Deg*(2^(Lev-1))+1:Deg*2^Lev;
epsilon = 1e-5;
eta = 5e-5;

xx = -1:0.01:1;
% begin of time advance
for t = 1:maxT
    % solve equation
    time = t*dt;
    Id = [];
    for i = 1:size(IHash,2)
        l1 = IHash{i};
        key = [l1(1),l1(2)];
        lr = l1(3);%Hash.(sprintf('i%g_',key));
        
        Id(Deg*(i-1)+[1:Deg]) = Deg*(lr-1)+[1:Deg];
    end
    MMeval = Meval(plot_start:num_In:end,Id);
    subplot(2,2,1)
%     plot(xx_node,MMeval*fval_MW,'r-o',xx_node,exactf(xx_node,time-dt),'b--','LineWidth',2)
    plot(xx_node,MMeval*fval_MW,'r-o',xx,exactf(xx,time-dt),'b--','LineWidth',2)
    title(['solution at time ',num2str(time-dt)])
%     axis([-1 1 -0.1 1])
    subplot(2,2,3)
    plotgrid(IHash,Lstart,Lend,0,'r-x');
    axis([-1 1 -1 1])
    title(['old grids at t= ',num2str(t)])
    
    pause(0.1)
    
    % pre-predict
    [Mat_tmp] = GlobalMatrix(Mat,IHash,Deg);
    % first solve it by Euler method for f^{p}
    fval_tmp = fval_MW + dt*( Mat_tmp*fval_MW );
    
    MarkedRefine = mark(fval_tmp,[1:size(IHash,2)*Deg],Deg,epsilon,'refine');
    
%     EndGrid = LeafSolIndex(end);
    
    % refinement
    [Hash,IHash,FineIndex] = RefineMesh(Hash,IHash,Deg,MarkedRefine,FineIndex);
%     size(IHash)
%     [size(IHash,2),FineIndex(end)]
%    FineIndex
%    Hash
    
    % construct new mat on the updated mesh
    [Mat_tmp] = GlobalMatrix(Mat,IHash,Deg);
    
    DoFs_new = size(Mat_tmp,1) - size(fval_tmp,1);
    fval_MW = [fval_MW;zeros(DoFs_new,1)];
    
    % compute f^{n+1}
    f1 = fval_MW + dt*( Mat_tmp*fval_MW );
    f2 = 3/4*fval_MW+1/4*f1+1/4*dt*( Mat_tmp*f1 );
    fval = 1/3*fval_MW+2/3*f2+2/3*dt*( Mat_tmp*f2 );
    
    fval_MW = fval;
    Id = [];
    for i = 1:size(IHash,2)
        l1 = IHash{i};
        key = [l1(1),l1(2)];
        lr = l1(3);%Hash.(sprintf('i%g_',key));
        
        Id(Deg*(i-1)+[1:Deg]) = Deg*(lr-1)+[1:Deg];
    end
    MMeval = Meval(plot_start:num_In:end,Id);
    
    subplot(2,2,2)
%     plot(xx_node,MMeval*fval_MW,'r-o',xx_node,exactf(xx_node,time),'b--','LineWidth',2);
    plot(xx_node,MMeval*fval_MW,'r-o',xx,exactf(xx,time),'b--','LineWidth',2);
    title(['solution at time ',num2str(time)])
%     axis([-1 1 -0.1 1])
% %     subplot(2,2,4)
% %     plotgrid(IHash,Lstart,Lend,0,'b-o');
% %     axis([-1 1 -1 1])
% %     title('new grids with refinement')
% %     pause(0.1)


    
    
% %     % coarsening
    checkIndex = Deg*(FineIndex-1)'+[1:Deg];
    checkIndex = sort(checkIndex(:));
    MarkedCoarsen = mark(fval_MW,checkIndex,Deg,eta,'coarse');
    [Hash,IHash,FineIndex,count_del] = CoarseMesh(Hash,IHash,Deg,MarkedCoarsen,FineIndex);
%      [size(IHash,2),FineIndex(end)]
% size(IHash)
% '=========='
    size(count_del)
%     FineIndex
    if size(count_del,2)>0
     index_del = Deg*(count_del-1)'+[1:Deg];
    index_del = index_del(:);
    fval_MW(index_del) = [];
    end
    
    subplot(2,2,4)
    plotgrid(IHash,Lstart,Lend,0,'b-o');
    axis([-1 1 -1 1])
    title('new grids with refinement and coarsen')
    pause(0.1)
    
%     if mod(t,20)==0
%         save(['IHash_t',num2str(t),'.mat'],'IHash','Deg','fval_MW','Hash')
%     end

end








