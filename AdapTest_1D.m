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

EpsMax = 1e-5;
EpsMin = EpsMax/10;

MaxIter = 10;
figure;
for iter = 1:MaxIter
    subplot(2,2,1);hold on
    idPlot = find(VecFlag>0);
    PlotGridsInd(Lstart,Lend,idPlot,VecFlag)
    title(['Grid 0 - ',num2str(size(idPlot,1))])
    % Input is f0val
    f1val = f0val;
    VecFlag1 = VecFlag;
    
    subplot(2,2,2);
    plot(xxplot,MM*f0val,'k-','LineWidth',2);hold on
    
    % checking for refinement
    Leaf4Refine = find(VecFlag > 1);
    Ind4Refine = Grid2Dof(Leaf4Refine,Deg);
    Val4Check = reshape(f0val(Ind4Refine),Deg,size(f0val(Ind4Refine),1)/Deg);
    tmp = MarkElem(Val4Check,'refine',EpsMax);
    IndGridRefine = Leaf4Refine(tmp);

    RefineLev = ceil(log2(IndGridRefine)) ;
    RefineCel = IndGridRefine - 1 - 2.^(RefineLev-1);
    
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
    PlotGridsInd(Lstart,Lend,idPlot,VecFlag1);
    title(['Grid p - ',num2str(size(idPlot,1))])
    hold off;
    
    % update coefficients
    AddDofInd = Grid2Dof(AddIndGrid,Deg);
    f1val(AddDofInd) = f_MW_full(AddDofInd);
    subplot(2,2,2);
    plot(xxplot,MM*f1val,'b-','LineWidth',2);
    
    % coarsen check
    Leaf4Coarse = find(VecFlag1 > 2);
    Ind4Coarse = Grid2Dof(Leaf4Coarse,Deg);
    Val4Check = reshape(f1val(Ind4Coarse),Deg,size(f1val(Ind4Coarse),1)/Deg);
    tmp = MarkElem(Val4Check,'coarse',EpsMin);
    IndGridCoarse = Leaf4Coarse(tmp);

    CoarseLev = ceil(log2(IndGridCoarse)) ;
    CoarseCel = IndGridCoarse - 1 - 2.^(CoarseLev-1);
    
    ixD = find(CoarseLev<=IniLev); % never delete lev less than initial
    CoarseLev(ixD) = [];
    CoarseCel(ixD) = [];
    
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
    f1val(DelIndDof) = 0;
    subplot(2,2,2);
    plot(xxplot,MM*f1val,'g--','LineWidth',2);
    hold off
    % update coefficients
    f0val = f1val; 
    VecFlag = VecFlag1;
    
    % plot grids
    subplot(2,2,4);hold on
    idPlot = find(VecFlag>0);
    PlotGridsInd(Lstart,Lend,idPlot,VecFlag);
    title(['Grid 1 - ',num2str(size(idPlot,1))])
    hold off;
end
% check the deepest layer
LeafLayer = find(VecFlag == 2);
IndMat = Grid2Dof(LeafLayer,Deg);

ix = find(abs(f0val(IndMat))>EpsMax);
ix = IndMat(ix);
% Place needs refinement
MarkedElem = unique(ceil(ix/Deg));
MarkedLev = ceil(log2(MarkedElem)) ;
MarkedCel = MarkedElem - 1 - 2.^(MarkedLev-1);

NewLev = MarkedLev+1;
NewCel = [2*MarkedCel,2*MarkedCel+1];

NewInd = (2.^(NewLev-1)+NewCel+1)';
NewInd = NewInd(:);

VecFlag(NewInd) = 2;
VecFlag(MarkedElem) = 1;

% plot the grid
IndNZ = find(VecFlag>0);
figure;hold on
for i = 1:size(IndNZ)
    if i == 1
        plot((Lend+Lstart)/2,0,'ko','MarkerFaceColor','k','MarkerSize',10)
    else
        tmp = IndNZ(i);
        LevT = ceil(log2(tmp));
        CelT = tmp-1-2^(LevT-1);

        hT = (Lend-Lstart)/2^(LevT-1);
        xT(i) = Lstart+CelT*hT+hT/2;
        if VecFlag(tmp) == 2
            plot(xT(i),-LevT,'bo','MarkerFaceColor','b','MarkerSize',10);
        else
            plot(xT(i),-LevT,'ro','MarkerFaceColor','r','MarkerSize',10);
        end
        
    end
    
end

% update fval
NewDofInd = Grid2Dof(NewInd,Deg);%(repmat(NewInd(:),1,Deg)-1)*Deg+[1:Deg];

f0val(NewDofInd) = f_MW_full(NewDofInd);

% coarsing
ixCoarse = find(VecFlag == 2);
ixDof = Grid2Dof(ixCoarse,Deg);
tmp = find(abs(f0val(ixDof))<1e-7);
ix = ixDof(tmp);
f0val(ix)
MarkedGrid = unique(ceil(ix/Deg));

MarkedLev = ceil(log2(MarkedGrid)) ;
MarkedCel = MarkedGrid - 1 - 2.^(MarkedLev-1);
MarkedLevP =  MarkedLev-1;

% repeat checking coefficients
% check the deepest layer
LeafLayer = find(VecFlag == 2);
IndMat = (repmat(LeafLayer,1,Deg)-1)*Deg+[1:Deg];
IndMat = sort(IndMat(:));
ix = find(abs(f0val(IndMat))>EpsMax);
ix = IndMat(ix);
% Place needs refinement
MarkedElem = unique(ceil(ix/Deg));
MarkedLev = ceil(log2(MarkedElem)) ;%+ 1;
MarkedCel = MarkedElem - 1 - 2.^(MarkedLev-1);

NewLev = MarkedLev+1;
NewCel = [2*MarkedCel,2*MarkedCel+1];

NewInd = (2.^(NewLev-1)+NewCel+1)';
NewInd = NewInd(:);

for i = 1:size(NewInd,1)
    tmp = NewInd(i);
    LevT = ceil(log2(tmp));
    CelT = tmp-1-2^(LevT-1);

    hT = (Lend-Lstart)/2^(LevT-1);
    xT = Lstart+CelT*hT+hT/2;
    plot(xT,-LevT,'gs','MarkerFaceColor','g','MarkerSize',10);
end
% check about following

Ind = 2^(LevT - 1) + CelT + 1;% only for LevT > 0, CelT = 0~2^(LevT-1)
LevT = ceil(log2(Ind));
CelT = Ind - 1 - 2^(LevT-1);

% Matrix

% Step 4. Time Advance
for t = 1:MaxT
    time = t*dt;
    
    fval_tmp = fval0 + dt*( Mat_tmp*fval_MW );
end

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
