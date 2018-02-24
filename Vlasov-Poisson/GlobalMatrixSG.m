function A_encode=GlobalMatrixSG(A,B,Hash)

% save the submatrix and apply A_encode to vector f
count=0;
%=====================================================
% The Index is grouped by the sum of Levels in V and X
% {Iv,Jv} together with Cell and Deg information will
%         determine the position of matrix A
% {Ix,Jx} together with Cell and Deg information will
%         determine the position of matrix B
% Note: SG needs Iv+Ix<=Lev and Jv+Jx<=Lev
%=====================================================
Deg = Hash.Deg;
Lev = Hash.Lev_1;
tic
for sum = 0:Lev % Loop for Total Sum
%     Count=((sum+1)^2*(sum+2)+sum*(sum+1)^2)2
    %---------------
    % Fix Ix+Iv=sum
    %---------------
    for Ix = 0:sum
        Iv = sum-Ix;
        
        for Jx = 0:sum
            for Jv = 0:sum-Jx
                
                [index_I,index_J,tmpA,tmpB]=ComputeGlobalIndex(Ix,Jx,Iv,Jv,Hash,A,B);    
%                 [Ix,Iv,Jx,Jv]
                % save matrices to files
                count=count+1;
                A_encode{count}.A1=tmpA;
                A_encode{count}.A2=tmpB;
                
                
                A_encode{count}.IndexI = index_I;
                A_encode{count}.IndexJ = index_J;
                
            end
        end
    end
%     disp('=====')
    %---------------
    % Fix Jx+Jv=sum
    %---------------
    for Jx=0:sum
        Jv=sum-Jx;
        
        % Avoid Duplicate Index
        for Ix=0:sum%-1
            for Iv=0:sum-Ix-1
%                 [Ix,Iv,Jx,Jv]
                [index_I,index_J,tmpA,tmpB]=ComputeGlobalIndex(Ix,Jx,Iv,Jv,Hash,A,B);
                
                % save matrices to files
                count=count+1;
                A_encode{count}.A1=tmpA;
                A_encode{count}.A2=tmpB;
                
                
                A_encode{count}.IndexI = index_I;
                A_encode{count}.IndexJ = index_J;
            end
        end
    end
end
% toc

% count
end

function [index_I,index_J,tmpA,tmpB]=ComputeGlobalIndex(Ix,Jx,Iv,Jv,Hash,A,B)

k = Hash.Deg;
lev = Hash.Lev_2;

index_I_x = IndexMapping(Ix,k);
index_J_x = IndexMapping(Jx,k);
index_I_v = IndexMapping(Iv,k);
index_J_v = IndexMapping(Jv,k);

tmpA = A(index_I_v,index_J_v);
tmpB = B(index_I_x,index_J_x);


Key_I = GenerateKey2D(index_I_v,index_I_x,k,lev);
Key_J = GenerateKey2D(index_J_v,index_J_x,k,lev);

index_I=zeros(size(Key_I,1),1);
for i=1:size(Key_I,1)
    index_I(i) = Hash.(sprintf('i%g_',Key_I(i,:)));
end

index_J=zeros(size(Key_J,1),1);
for i=1:size(Key_J,1)
    index_J(i) = Hash.(sprintf('i%g_',Key_J(i,:)));
end

end

function index=IndexMapping(I,Deg)
% This function is to find all the index from Lev (I-1) to Lev (I)
if I==0
    index=[1:Deg];
else
    index=Deg*2^(I-1)+1:Deg*2^(I);
end
end

function key=GenerateKey2D(Ix,Iy,Deg,Lev)
% First generate Key1dMesh
nx=[];px=[];kx=[];
for Lx=0:Lev
    for Px=0:2^max(0,Lx-1)-1
        for Kx=1:Deg
            nx=[nx;Lx];
            px=[px;Px];
            kx=[kx;Kx];
            
        end
    end
end
%
Key1dMesh=[nx,px,kx];

% Second generate the 2D Key
tmp_x=Key1dMesh(Ix,:);
tmp_y=Key1dMesh(Iy,:);
tmp_ix=reshape(repmat(tmp_x,1,size(tmp_y,1))',3,size(tmp_x,1)*size(tmp_y,1) )';
tmp_iy=repmat(tmp_y,size(tmp_x,1),1);
key=[tmp_ix(:,1),tmp_iy(:,1),tmp_ix(:,2),tmp_iy(:,2),tmp_ix(:,3)-1,tmp_iy(:,3)-1];

end
