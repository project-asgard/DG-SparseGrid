% plot solutions on Adaptive Grids
h=2^(-n);
xnode=[0:h:1];
xnode=[xnode(1:end-1);xnode(2:end)];
xnode=xnode(:);

Hash.print

I=1:length(xnode);
xxnode=xnode(I);
MMeval=Meval(I,:);

MapMW2DG=sparse(dof_sparse,size(xxnode,1)^2);

load('big map.dat');

All_index=allv(Hash);
All_index=sort(All_index);
Real_index=-ones(size(All_index,2),1);
Real_index(All_index+1,1)=[1:size(All_index,2)];

for i=1:length(big_map)
%     i
    nx=big_map(i,1);
    ny=big_map(i,2);
    px=big_map(i,3);
    py=big_map(i,4);
    kx=big_map(i,5);
    ky=big_map(i,6);
    
    key_1=[nx,px,kx];
    key_2=[ny,py,ky];

    M1=MMeval(:,Index1D(nx,px,kx+1,k));
    M2=MMeval(:,Index1D(ny,py,ky+1,k));

    tmp=M1*M2';
    MapMW2DG(Real_index(big_map(i,7)+1),:)=tmp(:);
    
end

disp('Done of MapMW2DG')

val=MapMW2DG'*sol_s;
val=reshape(val,size(xxnode,1),size(xxnode,1));
[xx,yy]=meshgrid(xxnode);
figure;
subplot(1,2,1)
mesh(xx,yy,val, 'FaceColor','interp',...
    'EdgeColor','interp');
% axis([0 1 0 1 -0.02 1 ])
subplot(1,2,2)
sol=exactu_2D(xx,yy);
mesh(xx,yy,exactu_2D(xx,yy), 'FaceColor','interp',...
    'EdgeColor','interp');

['The error for 2D approximation is:']
max(max(val-sol))


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
