% main code for LDG equation
%==================================================
%  df/dt + d/dx(x(1-x^2)f) = source
% A = (x(1-x^2)f,dw/dx) - <x(1-x^2)\hat{f},w>
% F{n+1}=F{n}+dt*A*F{n}+dt*S
% Flux is:: {funcCoef*u}+abs(funcCoef)*(1-alpha)[u]
% Note:
%   alpha = 0:: upwind  flux
%   alpha = 1:: central flux
% Parameters:
%   alpha:: flux choice as aobve
%   PlotType:	1--discontinous plot
%               2--continuous plot
%   Deg, Lev, CFL
%==================================================
clear all
close all
% clc

% Test
sigma = 0.1;
% f0 = @(x)( exp(-(x-0.36).^2/sigma^2) );
f0 = @(x)( exp(-(x-0.36).^2/sigma^2)+exp(-(x+0.36).^2/sigma^2) );
phi = @(x,t)( x.*exp(-t)./sqrt(1+(exp(-2*t)-1).*x.^2) );
exactf = @(x,t)(...
    (phi(x,t).*(1-phi(x,t).^2))./(x.*(1-x.^2)).*f0(phi(x,t)) ...
    );

f0 = @(x)(x.^5);
exactf = @(x,t)(f0(x));
funcCoef = @(x)(x.*(1-x.^2));

format short e
addpath(genpath(pwd))


Lev = 2;
Deg = 2;
num_plot = 3;
EndTime = 4;
PlotType = 1;
AdapTest = 0;

alpha = 0;

% Define domain = [Lstart,Lend] and Lmax = Lend-Lstart
Lstart = -1;
Lend = 1;
Lmax = Lend-Lstart;



%--Quadrature
quad_num=10;
%---------------

% compute the trace values
p_1 = legendre(-1,Deg);
p_2 = legendre(1,Deg);

[quad_x,quad_w]=lgwt(quad_num,-1,1);
p_val = legendre(quad_x,Deg);
Dp_val = dlegendre(quad_x,Deg);

%---------------------------
% Jacobi of variable x and v
% Define Matrices
%---------------------------
if AdapTest == 0
    tol_cel_num=2^(Lev);h=Lmax/tol_cel_num;
    dof_1D=Deg*tol_cel_num;
    
    CFL = 0.01;
    
    
    
    if Deg < 4
        dt = CFL*h^((Deg-1)/3)/2;
    else
        dt = CFL*h^((Deg)/3);
    end
    
    
    
    maxT = ceil(EndTime/dt)
    
    M_Adv = sparse(dof_1D,dof_1D);
    S = sparse(dof_1D,1);
    f0 = sparse(dof_1D,1);
    
    b = sparse(dof_1D,1);
    fexact = sparse(dof_1D,1);
    
    % generate 1D matrix for DG
    for L=0:tol_cel_num-1
        
        %---------------------------------------------
        % (funcCoef*q,d/dx p)
        %---------------------------------------------
        x0 = Lstart+L*h;
        x1 = x0+h;
        xi = quad_x*(x1-x0)/2+(x1+x0)/2;
        xmid= (x0+x1)/2;
        
        val=1/h*[Dp_val'*(quad_w.*funcCoef(xi).*p_val)];
        
        c = Deg*L+1:Deg*(L+1);
        
        M_Adv = M_Adv + sparse(c'*ones(1,Deg),ones(Deg,1)*c,val,dof_1D,dof_1D);
        
        
        val = sqrt(h)/2*[p_val'*(quad_w.*exactf(xi,0))];
        fexact(c)=fexact(c)+val;
        
        
        %----------------------------------------------
        % -<funcCoef*{q},p>
        %----------------------------------------------
        if (L >= 1) && (L < tol_cel_num-1)
            val_flux = [...
                p_1'*(funcCoef(x0)+abs(funcCoef(x0))*(1-alpha))*p_2;...
                p_1'*(funcCoef(x0)-abs(funcCoef(x0))*(1-alpha))*p_1;...
                -p_2'*(funcCoef(x1)+abs(funcCoef(x1))*(1-alpha))*p_2;...
                -p_2'*(funcCoef(x1)-abs(funcCoef(x1))*(1-alpha))*p_1;...
                ]/2/h;
            IndexV = [c'*ones(1,Deg);c'*ones(1,Deg);c'*ones(1,Deg);c'*ones(1,Deg)];
            IndexU = [ones(Deg,1)*(c-Deg);ones(Deg,1)*(c);ones(Deg,1)*(c);ones(Deg,1)*(c+Deg) ];
        elseif L == 0
            val_flux = [...
                p_1'*(funcCoef(x0)-abs(funcCoef(x0))*(1-alpha))*p_1;...
                -p_2'*(funcCoef(x1)+abs(funcCoef(x1))*(1-alpha))*p_2;...
                -p_2'*(funcCoef(x1)-abs(funcCoef(x1))*(1-alpha))*p_1;...
                ]*0.5/h;
            
            IndexV = [c'*ones(1,Deg);c'*ones(1,Deg);c'*ones(1,Deg)];
            IndexU = [ones(Deg,1)*(c);...
                ones(Deg,1)*(c);...
                ones(Deg,1)*(c+Deg) ];
            val_L = -p_1'*(funcCoef(x0)+abs(funcCoef(x0))*(1-alpha))*0;
        elseif L == tol_cel_num-1
            val_flux = [...
                p_1'*(funcCoef(x0)+abs(funcCoef(x0))*(1-alpha))*p_2;...
                p_1'*(funcCoef(x0)-abs(funcCoef(x0))*(1-alpha))*p_1;...
                -p_2'*(funcCoef(x1)+abs(funcCoef(x1))*(1-alpha))*p_2;...
                ]*0.5/h;
            
            IndexV = [c'*ones(1,Deg);c'*ones(1,Deg);c'*ones(1,Deg)];
            IndexU = [ones(Deg,1)*(c-Deg); ...
                ones(Deg,1)*(c);...
                ones(Deg,1)*(c)];
            val_R = p_2'*(funcCoef(x1)-abs(funcCoef(x1))*(1-alpha))*0;
        end
        M_Adv = M_Adv + sparse(IndexV,IndexU,val_flux,dof_1D,dof_1D);
        
        val = sqrt(h)/2*[p_val'*(quad_w.*exactf(xi,0))];
        f0(c) = val;
        
        
    end
    
    
    
    
    [quad_x,quad_w]=lgwt(num_plot,-1,1);
    
    p_val = legendre(quad_x,Deg);
    for L=0:tol_cel_num-1
        %---------------------------------------------
        % Generate the coefficients for DG bases
        %---------------------------------------------
        
        Iu = [Deg*L+1:Deg*(L+1)];
        
        Iv = [num_plot*L+1:num_plot*(L+1)];
        
        x0 = Lstart+L*h;
        x1 = x0+h;
        
        if PlotType == 1
            xi = [x0,quad_x(2:end-1)'*(x1-x0)/2+(x1+x0)/2,x1];
        elseif PlotType == 0
            xi = quad_x*(x1-x0)/2+(x1+x0)/2;
        end
        
        
        Meval(Iv,Iu)=sqrt(1/h)*p_val;
        x_node(Iv,1)=xi;
        
    end
    
    Mat = M_Adv;
    
    b = S;
    
    % convect matrix to MWDG
    FMWT = OperatorTwoScale(Deg,2^Lev);
    %         Mat = FMWT*Mat*FMWT';
    %         b = FMWT*b;
    Meval = Meval*FMWT';
    %
    f0 = FMWT*f0;
    for i = 0:Lev
        if i == 0
            startP = 1;
        else
            startP = Deg*(2^max(i-1,0))+1;
        end
        endP = Deg*2^i;
        fcell{i+1} = f0(startP:endP);
        
    end
    %
    %         save(['Mat_Lev',num2str(Lev),'_Deg',num2str(Deg),'.mat']);
    
    % checked of projection
    plot(x_node,Meval*f0,'r-o',x_node,exactf(x_node,0),'b--','LineWidth',2)
    legend({'solution f','flux x(1-x^2)f'})
    title(['time at ',num2str(0)])
    
    return
    
    
    
    [quad_x,quad_w]=lgwt(num_plot,-1,1);
    total_particle = 0;
    ffval = Meval*f0;
    L2_stability = 0;
    for i = 1:num_plot
        total_particle =  total_particle+quad_w(i)*h/2*sum(ffval(i:num_plot:end));
        L2_stability = L2_stability+quad_w(i)*h/2*sum(ffval(i:num_plot:end).^2);
    end
    [total_particle L2_stability]
    
    tp(1) = total_particle;
    Lp(1) = L2_stability;
    
    figure
    for t = 1:maxT
        %     t
        time = t*dt;
        
        % if source is time-independent
        f1 = f0 + dt*( Mat*f0+b );
        f2 = 3/4*f0+1/4*f1+1/4*dt*(Mat*f1+b);
        fval = 1/3*f0+2/3*f2+2/3*dt*(Mat*f2+b);
        % if source is time-dependent
        %     f1 = f0 + dt*( Mat*f0+b*exp(time-dt) );
        %     f2 = 3/4*f0+1/4*f1+1/4*dt*(Mat*f1+b*exp(time));
        %     fval = 1/3*f0+2/3*f2+2/3*dt*(Mat*f2+b*exp(time-dt/2));
        
        f0 = fval;
        
        plot(x_node,Meval*f0,'r-o',x_node,x_node.*(1-x_node.^2).*Meval*(f0),'b-<')
        legend({'solution f','flux x(1-x^2)f'})
        title(['time at ',num2str(time)])
        pause (0.1)
        
        total_particle = 0;
        L2_stability = 0;
        ffval = Meval*f0;
        for i = 1:num_plot
            
            total_particle =  total_particle+...
                quad_w(i)*h/2*sum(ffval(i:num_plot:end));
            L2_stability = L2_stability+quad_w(i)*h/2*sum(ffval(i:num_plot:end).^2);
        end
        tp(t+1) = total_particle;
        Lp(t+1) = L2_stability;
        
        % saving data
        if abs(time-0.5)<=dt || abs(time-1)<=dt || abs(time-1.5)<=dt ...
                || abs(time-2)<=dt || abs(time-2.5)<=dt || abs(time-3)<=dt ...
                || abs(time-3.5)<=dt || abs(time-4)<=dt || abs(time-4.5)<=dt ...
                || abs(time-5)<=dt || abs(time-5.5)<=dt || abs(time-6)<=dt
            
            save(['Damp_ALPHA',num2str(alpha),'_Deg',num2str(Deg),'_Lev',num2str(Lev),'_End',num2str(time),'.mat'])
        end
    end
    
    val = Meval*f0-exactf(x_node,time);
    
    
    fL2 = 0; fLinf = max(abs(val));
    
    ffval = Meval*f0;
    for i = 1:num_plot
        fL2 = fL2 + quad_w(i)*h/2*sum(val(i:num_plot:end).^2);
    end
    [sqrt(fL2) fLinf]
    
    figure;
    plot(x_node,Meval*f0,'r-o',x_node,exactf(x_node,time),'r--','LineWidth',2);
    
    figure;
    plot(tp,'r-o'); hold on; plot(Lp,'b-o')
    
    
    
elseif AdapTest == 1
    
    
    CFL = 0.01;
    if Deg < 4
        dt = CFL*(Lmax/2^Lev)^((Deg-1)/3)/2;
    else
        dt = CFL*(Lmax/2^Lev)^((Deg)/3);
    end
    maxT = ceil(EndTime/dt)
    
    load(['Mat_Lev',num2str(10),'_Deg',num2str(Deg),'.mat'],'Mat','f0','Meval','x_node');
    
    [Hash,IHash,FineIndex] = HashTable1D(Lev);
    DoFs = Deg*2^Lev;
    
    Mat_tmp = Mat(1:DoFs,1:DoFs);
    
    
    for i = 1:size(IHash,2)
        l1 = IHash{i};
        Id(Deg*(i-1)+[1:Deg]) = Deg*(l1(3)-1)+[1:Deg];
    end
    
    num_In = 2^5;%2^(10-Lev);
    plot_start = 2^(10-Lev);
    MMeval = Meval(plot_start:num_In:end,Id);
    xx_node = x_node(plot_start:num_In:end);
    fval_MW = f0(1:DoFs);
    figure;
    plot(xx_node,MMeval*fval_MW,'r-o',xx_node,exactf(xx_node,0),'b--','LineWidth',2)
    legend({'solution f','flux x(1-x^2)f'})
    title(['time at ',num2str(0)])
    
    %
    LeafSolIndex = Deg*(2^(Lev-1))+1:Deg*2^Lev;
    epsilon = 1e-6;
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
        plot(xx_node,MMeval*fval_MW,'r-o',xx,exactf(xx,time-dt),'b--','LineWidth',2);%,
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
        %         fval_MW([(FineIndex-1)*Deg+1;(FineIndex-1)*Deg+2;(FineIndex-1)*Deg+3;FineIndex*Deg])
        
        subplot(2,2,4)
        plotgrid(IHash,Lstart,Lend,0,'b-o');
        axis([-1 1 -1 1])
        title('new grids with refinement and coarsen')
        pause(0.1)
        
        if mod(t,20)==0
            save(['Damp_IHash_t',num2str(t),'.mat'],'IHash','Deg','fval_MW','Hash')
        end
        
    end
    
    
end



