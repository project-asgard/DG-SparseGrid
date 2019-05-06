% This is 1D adaptivity test for hyperbolic equation
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


% First set MaxLev = 12
MaxLev = 6;
Num_MaxNode = 2^MaxLev;
MaxDof = Deg*Num_MaxNode;

IniLev = 4;
Num_IniNode = 2^IniLev;
IniDof = Deg*Num_IniNode;

CFL = 0.1;
h = (Lend-Lstart)/2^MaxLev;
dt = .1*h;%^((Deg)/3);
maxT = ceil(EndTime/dt);

% Step 1. Calculate a large grad operator
Mat = MatrixGrad(MaxLev,Deg,Lstart,Lend,1,funcCoef );


% Step 2. List of FLAG
VecFlag = sparse(Num_MaxNode,1);
VecFlag(1:Num_IniNode) = 1;
VecFlag(Num_IniNode - 2^(IniLev-1)+1:Num_IniNode) = 3; % This denotes the deepest layer


% Step 3. Initial Condition
quad_num = Deg;
[quad_x,quad_w] = lgwt(quad_num,-1,1);
pV = lin_legendre2(quad_x,Deg)*1/sqrt(h);


fval = sparse(MaxDof,1);
fval0 = sparse(MaxDof,1);
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

fval0(1:IniDof,1) = f_MW_full(1:IniDof);
% Mat = FMWT* Mat*FMWT';

load('Mat_Lev10_Deg2_v2.mat','Mat')
IdOn = find(VecFlag>0);
IdDofOn = Grid2Dof(IdOn,Deg);
% A_encode{1}.A = Mat(IdDofOn,IdDofOn);
% A_encode{1}.IndexI = IdDofOn;
% A_encode{1}.IndexJ = IdDofOn;

Iu = IdDofOn*ones(1,size(IdDofOn,1));
Aij = Mat(IdDofOn,IdDofOn);
Mat_A = sparse(Iu,Iu',Aij);

EpsMax = 1e-5;
EpsMin = EpsMax/100;
% begin of time loop
count = 1;

fval_tmp = sparse(MaxDof,1);
fval_tmp0 = sparse(MaxDof,1);
figure;
time = 0;
for iter = 1:1000%MaxIter
    iter
    time = time + dt;
    VecFlag1 = VecFlag;
%     B_encode = A_encode;
    
    subplot(1,2,1);
    plot(xxplot,MM*fval0,'b-','LineWidth',2);hold on;
    idPlot = find(VecFlag>0);
    PlotGridsInd(Lstart,Lend,idPlot,VecFlag,2);
    
    % pre-computing f1
    %     fval_tmp = fval0 + dt*( ApplyA(B_encode,fval0,VecFlag) );
    
    val_tmp = Mat_A*fval0(IdDofOn);
    fval_tmp0(IdDofOn) = fval0(IdDofOn) + dt*val_tmp;
    
    %     subplot(1,3,2);
    plot(xxplot,MM*fval_tmp0,'r--','LineWidth',2);
    title('Input Grids')
    hold off;
    %
    % checking for refinement
    Leaf4Refine = find(VecFlag > 1);
    Ind4Refine = Grid2Dof(Leaf4Refine,Deg);
    Val4Check = reshape(fval_tmp0(Ind4Refine),Deg,size(fval_tmp(Ind4Refine),1)/Deg);
    tmp = MarkElem(Val4Check,'refine',EpsMax);
    IndGridRefine = Leaf4Refine(tmp);
    
    RefineLev = ceil(log2(IndGridRefine)) ;
    RefineCel = IndGridRefine - 1 - 2.^(RefineLev-1);
    
    ixd = find(RefineLev>=MaxLev);
    RefineLev(ixd) = [];
    RefineCel(ixd) = [];
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
    
    % update coefficients
    AddDofInd = Grid2Dof(AddIndGrid,Deg);
    fval0_tmp = fval0;
    
    % update matrix
    IdOn = find(VecFlag1>0);
    
    IdDofOn = Grid2Dof(IdOn,Deg);
    Mat_B = Mat(IdDofOn,IdDofOn);
    f1 = fval0(IdDofOn)+dt*Mat_B*fval0(IdDofOn);
    f2 = 3/4*fval0(IdDofOn)+1/4*f1+1/4*dt*(Mat_B*f1 );
    val_tmp = 1/3*fval0(IdDofOn)+2/3*f2+2/3*dt*( Mat_B*f2 );
    fval_tmp(IdDofOn) = val_tmp;
%     fval_tmp(IdDofOn) = fval0(IdDofOn) + dt*val_tmp;
    
    
    subplot(1,2,2);
    plot(xxplot,MM*fval_tmp,'b-','LineWidth',2);
    hold on;
    idPlot = find(VecFlag1>0);
    PlotGridsInd(Lstart,Lend,idPlot,VecFlag1,2);
    title(['Pre-Refinement ',num2str( size(idPlot,1))])
    %     hold off;
    
    % coarsen check
    Leaf4Coarse = find(VecFlag1 > 2);
    Ind4Coarse = Grid2Dof(Leaf4Coarse,Deg);
    Val4Check = reshape(fval_tmp(Ind4Coarse),Deg,size(fval_tmp(Ind4Coarse),1)/Deg);
    tmp = MarkElem(Val4Check,'coarse',EpsMin);
    IndGridCoarse = Leaf4Coarse(tmp);
    
    CoarseLev = ceil(log2(IndGridCoarse)) ;
    CoarseCel = IndGridCoarse - 1 - 2.^(CoarseLev-1);
    
%     ixD = find(CoarseLev<=IniLev); % never delete lev less than initial
%     CoarseLev(ixD) = [];
%     CoarseCel(ixD) = [];
    
    DelIndGrid = (2.^(CoarseLev-1)+CoarseCel+1)';
    DelIndGrid = DelIndGrid(:);
    
    % delete grid
    VecFlag1(DelIndGrid) = 0;
    
    % update VecFlag
    LevUp = CoarseLev-1;
    CelUp = ceil((CoarseCel-1)/2);
    IndUp = (2.^(LevUp-1)+CelUp+1)';
    VecFlag1(IndUp) = 2;
    
    DelIndDof = Grid2Dof(DelIndGrid,Deg);
    fval_tmp(DelIndDof) = 0;
    % update the deepest grids
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
    
    %     subplot(1,3,2);
    plot(xxplot,MM*fval_tmp,'r--','LineWidth',2);
    hold on;
    plot(xxplot,exactf(xxplot,time),'g-','LineWidth',2)
    idPlot = find(VecFlag1>0);
    PlotGridsInd(Lstart,Lend,idPlot,VecFlag1,2);
    title(['Post-Coarsen ',num2str( size(idPlot,1))])
    hold off;
    %    [size(A_encode), size(B_encode)]
    pause(0.1);
    
    % update grids and coefficients
    Mat_A = Mat_B;
    fval0 = fval_tmp;
    VecFlag = VecFlag1;
end

function ff = ApplyA(A,f,VecFlag)
nz = size(A,2);
Dof = size(f,1);
ff = sparse(Dof,1);
for i = 1:nz
    indexI = A{i}.IndexI;
    indexJ = A{i}.IndexJ;
    f_tmp = f(indexJ);
    ix = find(VecFlag(indexJ)==0);
    f_tmp = 0;
    tmp = A{i}.A*f(indexJ);
    ix = find(VecFlag(indexI)==0);
    tmp( (ix) ) = 0;
    ff(indexI) = ff(indexI) + tmp;%A{i}.A*f(indexJ);
    clear tmp
end
end