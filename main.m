% This code is designed to solve the Poisson problem
% by Hp-Weak Galerkin Methods
clc
close all
clear all
format short e

isplotmesh=0;

%---Give the value for the following parameters----------------------------
testproblem=1;
polydeg = 2;
%--------------------------------------------------------------------------
iter=1;
a=0;b=1;


for nn=1:1
    n=1%2^nn;
    
%% Generate the uniform triangular mesh
blog('Generate the uniform triangles');
[meshData] = meshGenerateSlantShaped(a, b, a, b, n, n);
elog;
%-------Plot mesh----------------------------------------------------------
if isplotmesh==1
triplot(meshData.T',meshData.P(1,:),meshData.P(2,:),'LineWidth',1)
end

%% Solve the second order elliptic equation
uhp=hpWGsolver(meshData,polydeg);


[h,cond,x,Querror,QuL2,uL2inf,ExactU,uerror,uL2,uH2error]=...
    assembleStiffMatrix_MWG_PoissonV1(meshData,testproblem);

% [h,cond,x,Querror,QuL2,uL2inf,ExactU,uerror,uL2,uH2error]=...
%     assembleStiffMatrix_MWG_PoissonV2(meshData,testproblem);

QH1error(iter)=Querror;
QL2error(iter)=QuL2;
L2inf(iter)=uL2inf;
H1error(iter)=uerror;
L2error(iter)=uL2;
H2error(iter)=uH2error;

CC(iter)=cond

h(iter)=E;

% [h',QH1error',QL2error',L2inf',H1error',L2error',H2error']
[h',QH1error',QL2error']
iter=iter+1

end
