function Con2D=Connect2D(Lev,HASH,HASHInv,gridType)
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

% ----------------------
% option to sort index_J
% ----------------------
do_sort = 1;

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


Lev_tmp = Lev;
if gridType == 'FG'
        Lev_tmp = Lev + Lev;
elseif gridType == 'SG'
        Lev_tmp = Lev;
end

for ii=1:nHash
    
    ll=HASHInv{ii};
    
    % 1D indices by the HashInv
    I1=ll(5);I2=ll(6);
    
    % find the connectivity from 1D connected mesh
    [i,j,val]=find(Con1D(I1,:));
    
    % Get (m1,cell1) from Key1DMesh
    LevCell1=Key1dMesh(j,:);
    
    [i,j,val]=find(Con1D(I2,:));
    
    % Get (m2,cell2) from Key1DMesh
    LevCell2=Key1dMesh(j,:);

    
    index_J=[];
    for i1=1:size(LevCell1,1)
        for i2=1:size(LevCell2,1)
            
                    
            if LevCell1(i1,1)+LevCell2(i2,1)<=Lev_tmp%Lev % check whether m1+m2<=Lev % modified by Lin
                
                key=[LevCell1(i1,1) LevCell2(i2,1) LevCell1(i1,2) LevCell2(i2,2)];
                index_J = [index_J, HASH.(sprintf(hash_format,key))];
                
            end
            
        end
    end

    if (do_sort),
            index_J = sort(index_J);
    end;
    
    Con2D{ii} = index_J;
    
end

%% Plotting for validation
%----------------------------
% matrix for 2D connectivity
%----------------------------
% Con2D_full=sparse(size(HASHInv,2));
% for i=1:size(Con2D,2)
%     Con2D_full(i,Con2D{i})=1;
% end

% figure;spy(Con2D_full)
% return
end

