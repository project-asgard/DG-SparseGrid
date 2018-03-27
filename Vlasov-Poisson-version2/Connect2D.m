function Con2D=Connect2D(Lev,HASH,HASHInv)
%================================================================
% This code is to generate the 2D connectivity based on Lev,HASH
% Here, we consider the maximum connectivity
% which including all overlapping cells, neighbor cells, and
% the periodic boundary cells
% Input: Lev, HASH, HASHInv
% Output: Con2D
%================================================================
%% Step 1. generate 1D connectivity
%=============================================================
% construct 1D connectivity
% output: the nonzero parts for matrices
% we consider the most numbers for connected case
% (ignoring Deg)
%=============================================================
global hash_format
Con=Connect1D(Lev);

%% Step 2. All possible combinations for 1D mesh
%--------------------------------------
% all possible combinations for 1D mesh
% output: 1D index--> (Lev,Cell)
% (ignoring Deg)
%--------------------------------------
nx=[];px=[];
for Lx=0:Lev
    for Px=0:2^max(0,Lx-1)-1
            nx=[nx;Lx];
            px=[px;Px];
    end
end
Key1dMesh=[nx,px];

%% Step 3. 2D connectivity
%--------------------------------------
% construct 2D connectivity 
% (ignoring Deg)
%--------------------------------------
for ii=1:HASH.dof
    ll=HASHInv{ii};
    
    n1=ll(1);p1=ll(3);
    n2=ll(2);p2=ll(4);
    
    I1=LevCell2index(n1,p1);
    J1=Con(I1,:);
    I2=LevCell2index(n2,p2);
    J2=Con(I2,:);
    
    [i,j,val]=find(Con(I1,:));
    LevCell1=Key1dMesh(j,:);
    
    [i,j,val]=find(Con(I2,:));
    LevCell2=Key1dMesh(j,:);
 
    index_J=[];
    for i1=1:size(LevCell1,1)
        for i2=1:size(LevCell2,1)
            if LevCell1(i1,1)+LevCell2(i2,1)<=Lev
                key=[LevCell1(i1,1) LevCell2(i2,1) LevCell1(i1,2) LevCell2(i2,2)];
%                 index_J = [index_J, HASH.(sprintf('i%g_',key))];
                    index_J = [index_J, HASH.(sprintf(hash_format,key))]; % suggested by Ed
            end
        end
    end
    
Con2D{ii}=index_J;

end

%% Plotting for validation
%----------------------------
% matrix for 2D connectivity
%----------------------------
% Con2D_full=sparse(HASH.dof);
% for i=1:size(Con2D,2)
%     Con2D_full(i,Con2D{i})=1;
% end
% 
% figure;spy(Con2D_full)
% return
end

function Con=Connect1D(Lev)
%=============================================================
% construct 1D connectivity
% output: the nonzero parts for matrices
% we consider the most numbers for connected case
% (ignoring Deg)
%=============================================================
Con=sparse(2^Lev,2^Lev);
for Lx=0:Lev
    for Px=0:2^max(0,Lx-1)-1
        
        I=LevCell2index(Lx,Px);
        J=LevCell2index(Lx,[max(Px-1,0):min(Px+1,2^max(0,Lx-1)-1)]);
        % diagnal is connected
        Con(I,J)=1;
        Con(J,I)=1;
        
        % periodic boundary is connected
        if Px==0
            tmp_end=LevCell2index(Lx,2^max(0,Lx-1)-1);
            Con(I,tmp_end)=1;
            Con(tmp_end,I)=1;
        elseif Px==2^max(0,Lx-1)-1
            tmp_begin=LevCell2index(Lx,0);
            Con(I,tmp_begin)=1;
            Con(tmp_begin,I)=1;
        end
        
        for Ly=Lx+1:Lev
            nL=Ly-Lx;
            
            % the overlapping part, one cell Left, and one cell Right
            if Lx>0
                Py=[max(2^(nL)*Px-1,0):min(2^(nL)*Px+2^(nL),2^(Ly-1)-1)];
            elseif Lx==0
                Py=[max(2^(nL-1)*Px-1,0):min(2^(nL-1)*Px+2^(nL-1),2^(Ly-1)-1)];
            end
            
            % periodic boundary is connected
            Py=[0,Py,2^max(0,Ly-1)-1];
            
            J=LevCell2index(Ly,Py);
            
            Con(I,J)=1;
            Con(J,I)=1;
        end
        
    end
end

end






