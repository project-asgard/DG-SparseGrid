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

Con1D = Connect1D(Lev);
nHash = numel(HASHInv);

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

Con2D = {}; % Empty cell array.

for ii=1:nHash
    
    ll=HASHInv{ii};
    
    % 1D indices from the HashInv
    I1=ll(3,1);I2=ll(3,2);
    
    % find the connected information
    J1=Con1D(I1,:);
    J2=Con1D(I2,:);
    
    % find the connectivity from 1D connected mesh
    [i,j,val]=find(Con1D(I1,:));
    
    % Get (m1,cell1) from Key1DMesh
    LevCell1=Key1dMesh(j,:);
    
    [i,j,val]=find(Con1D(I2,:));
    
    % Get (m2,cell2) from Key1DMesh
    LevCell2=Key1dMesh(j,:);
    
    index_J=[];
    key=[];
    for i1=1:size(LevCell1,1)
        
        lev1 = LevCell1(i1,1);
        cell1 = LevCell1(i1,2);
        
        for i2=1:size(LevCell2,1)
            
            lev2 = LevCell2(i2,1);
            cell2 = LevCell2(i2,2);
            
            if lev1+lev2 <= Lev % check whether m1+m2<=Lev
                
                key(:,1) = [lev1,cell1]; % dim1 (lev,cell)
                key(:,2) = [lev2,cell2]; % dim2 (lev,cell)
                keystr = getHashKeyStr(key);
                
                index_J = [index_J, HASH.(keystr)];
                
            end
            
        end
    end
    
    Con2D{ii} = index_J;
    
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

