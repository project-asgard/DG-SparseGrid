% plot solutions on Adaptive Grids
h=2^(-n);
xnode=[0:h:1];
xnode=[xnode(1:end-1);xnode(2:end)];
xnode=xnode(:);

Hash.print
I=1:length(xnode);
% I=1:floor(length(xnode)/20):length(xnode);
xxnode=xnode(I);
MMeval=Meval(I,:);

MapMW2DG=sparse(dof_sparse,size(xxnode,1)^4);

load('big map.dat');

All_index=allv(Hash);
All_index=sort(All_index);
Real_index=-ones(size(All_index,2),1);
Real_index(All_index+1,1)=[1:size(All_index,2)];

for i=1:length(big_map)
%     i
    n1=big_map(i,1);
    n2=big_map(i,2);
    n3=big_map(i,3);
    n4=big_map(i,4);
    
    i1=big_map(i,5);
    i2=big_map(i,6);
    i3=big_map(i,7);
    i4=big_map(i,8);    
    
    k1=big_map(i,9);
    k2=big_map(i,10);
    k3=big_map(i,11);
    k4=big_map(i,12);
    
    M1=MMeval(:,Index1D(n1,i1,k1+1,k));
    M2=MMeval(:,Index1D(n2,i2,k2+1,k));
    M3=MMeval(:,Index1D(n3,i3,k3+1,k));
    M4=MMeval(:,Index1D(n4,i4,k4+1,k));

%     tmp=kron(M1,kron(M2,kron(M3,M4)));
   tmp=kron(kron(kron(M1,M2),M3),M4);

    MapMW2DG(Real_index(big_map(i,13)+1),:)=tmp(:);
    
end

disp('Done of MapMW2DG')

val=MapMW2DG'*sol_s;
[x1,x2,x3,x4]=ndgrid(xxnode);

sol=exactu_4D(x1(:),x2(:),x3(:),x4(:));

max(val-sol(:))
% max(max(max(max(val-sol))))

function z=Index1D(level,cell,degree,k)
% Determine the 1D index
% [level,cell,degree]
if level==0
    level=0;
else
    level=2^(level-1);
end
z=k*level+k*cell+degree;
% z
end
