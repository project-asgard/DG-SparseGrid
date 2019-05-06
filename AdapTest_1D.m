% This is a 1D test for adaptivity test
close all
clear
clc

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
Deg = 2;
num_plot = Deg;

Lstart = -1;
Lend = 1;
EndTime = 3;

% Step 1. Calculate a large grad operator
% First set MaxLev = 12
MaxLev = 10;
Num_MaxNode = 2^MaxLev;
MaxDof = Deg*Num_MaxNode;

IniLev = 3;
Num_IniNode = 2^IniLev;
IniDof = Deg*Num_IniNode;

CFL = 0.01;
h = (Lend-Lstart)/2^MaxLev;
dt = CFL*h^((Deg)/3);
maxT = ceil(EndTime/dt);

% Step 2. List of FLAG
VecFlag = sparse(Num_MaxNode,1);
VecFlag(1:Num_IniNode) = 1;
VecFlag(Num_IniNode - 2^(IniLev-1)+1:Num_IniNode) = 3; % This denotes the deepest layer


% Step 3. Initial Condition
quad_num = Deg;
[quad_x,quad_w] = lgwt(quad_num,-1,1);
pV = lin_legendre2(quad_x,Deg)*1/sqrt(h);


fval = sparse(MaxDof,1);
f0val = sparse(MaxDof,1);
for Num_RealSpaceCell=0:2^MaxLev-1
    
    x0 = Lstart+Num_RealSpaceCell*(Lend-Lstart)/2^MaxLev;
    x1 = x0+(Lend-Lstart)/2^MaxLev;
    xi = quad_x*(x1-x0)/2+(x1+x0)/2;
    
    val = (h)/2*[pV'*(quad_w.*exactf(xi,0))];
    
    c = Deg*Num_RealSpaceCell+[1:Deg];
    fval(c,1) = val;
    
    xxplot(c,1) = xi;
    MMeval(c,c) = pV;
    
end

FMWT = OperatorTwoScale(Deg,2^MaxLev);
f_MW_full =  FMWT*fval;
MM = MMeval*FMWT';

f0val(1:IniDof,1) = f_MW_full(1:IniDof);
% % plot(xxplot,MM*f_MW_full)
% % hold on;
% % plot(xxplot,MM*f0val)

EpsMax = 1e-10;
EpsMin = EpsMax/10;

MaxIter = 20;
figure;

for iter = 1:MaxIter
    
    subplot(2,2,1);hold on
    idPlot = find(VecFlag>0);
    PlotGridsInd(Lstart,Lend,idPlot,VecFlag,1)
    title(['Grid 0 - ',num2str(size(idPlot,1)),' iter ',num2str(iter)])
    axis([-1 1 -10 0])
    hold off;
    
    % Input is f0val
    f1val = f0val;
    VecFlag1 = VecFlag;
    
    subplot(2,2,2);
    plot(xxplot,MM*f0val,'k-','LineWidth',2);
    hold on
    
    % checking for refinement
    Leaf4Refine = find(VecFlag1 > 1);
    Ind4Refine = Grid2Dof(Leaf4Refine,Deg);
    Val4Check = reshape(f0val(Ind4Refine),Deg,size(f0val(Ind4Refine),1)/Deg);
    tmp = MarkElem(Val4Check,'refine',EpsMax);
    IndGridRefine = Leaf4Refine(tmp);
    
    RefineLev = ceil(log2(IndGridRefine)) ;
    RefineCel = IndGridRefine - 1 - 2.^(RefineLev-1);
    
    ix_tmp = find(RefineLev>=MaxLev);
    RefineLev(ix_tmp) = [];
    RefineCel(ix_tmp) = [];
    
    % temporary adding Lev and Cel
    AddLev = RefineLev+1;
    AddCel = [2*RefineCel,2*RefineCel+1];
    
    AddIndGrid = (2.^(AddLev-1)+AddCel+1)';
    AddIndGrid = AddIndGrid(:);
    
    idNotNeed = find(VecFlag1(AddIndGrid) >0);
    AddIndGrid(idNotNeed) = [];
    
    % update grids and flags
    VecFlag1(AddIndGrid) = 3;
    VecFlag1(IndGridRefine) = 1;
    
    subplot(2,2,3);hold on
    idPlot = find(VecFlag1>0);
    PlotGridsInd(Lstart,Lend,idPlot,VecFlag1,1);
    title(['Grid p - ',num2str(size(idPlot,1))])
    axis([-1 1 -10 0])
    hold off;
    
    % update coefficients
    AddDofInd = Grid2Dof(AddIndGrid,Deg);
    f1val(AddDofInd) = f_MW_full(AddDofInd);
    subplot(2,2,2);
    plot(xxplot,MM*f1val,'b-','LineWidth',2);
%     hold off;
    
    % coarsen check
    Leaf4Coarse = find(VecFlag1 > 2);
    Ind4Coarse = Grid2Dof(Leaf4Coarse,Deg);
    Val4Check = reshape(f1val(Ind4Coarse),Deg,size(f1val(Ind4Coarse),1)/Deg);
    tmp = MarkElem(Val4Check,'coarse',EpsMin);
    IndGridCoarse = Leaf4Coarse(tmp);
    
    CoarseLev = ceil(log2(IndGridCoarse)) ;
    CoarseCel = IndGridCoarse - 1 - 2.^(CoarseLev-1);
    
    % %     ixD = find(CoarseLev<=IniLev); % never delete lev less than initial
    % %     CoarseLev(ixD) = [];
    % %     CoarseCel(ixD) = [];
    
    DelIndGrid = (2.^(CoarseLev-1)+CoarseCel+1)';
    DelIndGrid = DelIndGrid(:);
    
    % delete grid
    VecFlag1(DelIndGrid) = 0;
    
    % update VecFlag
    LevUp = CoarseLev-1;
    CelUp = ceil((CoarseCel-1)/2);
    IndUp = (2.^(LevUp-1)+CelUp+1)';
    VecFlag1(IndUp) = 2;
    
    % check whether no chidren
    for Num_Coarse = 1:size(LevUp,1)
        LevDown = LevUp(Num_Coarse)+1;
        CelDown = [2*CelUp(Num_Coarse),2*CelUp(Num_Coarse)+1];
        CheckIndGrid = (2.^(LevDown-1)+CelDown+1)';
        if VecFlag1(CheckIndGrid(1)) == 0 && VecFlag1(CheckIndGrid(2)) == 0
            LevUp2 = LevUp(Num_Coarse);
            CelUp2 = CelUp(Num_Coarse);
            IndUp2 = (2.^(LevUp2-1)+CelUp2+1)';
            VecFlag1(IndUp2) = 3;
            %         1111
        end
    end
    
    
    DelIndDof = Grid2Dof(DelIndGrid,Deg);
    f1val(DelIndDof) = 0;
    subplot(2,2,2);
    plot(xxplot,MM*f1val,'g--','LineWidth',2);
    hold off
    
    % plot grids
    subplot(2,2,4);hold on
    idPlot = find(VecFlag1>0);
    PlotGridsInd(Lstart,Lend,idPlot,VecFlag1,1);
    title(['Grid 1 - ',num2str(size(idPlot,1))])
    hold off;
    axis([-1 1 -10 0])
    if norm(VecFlag-VecFlag1) == 0
        break
    end
    
    
    %     % plot grids
    %     subplot(2,2,4);hold on
    %     idPlot = find(VecFlag>0);
    %     PlotGridsInd(Lstart,Lend,idPlot,VecFlag);
    %     title(['Grid 1 - ',num2str(size(idPlot,1))])
    %     hold off;
    
    pause(0.3)
    
    % update coefficients
    f0val = f1val;
    VecFlag = VecFlag1;
end

iter
