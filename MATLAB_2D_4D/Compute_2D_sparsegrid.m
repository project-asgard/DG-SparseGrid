% Computing on the sparse grid for Poisson Eq
% 2-dimensional calculation
%------------------------------------------------
% A_s: Matrix
% b_s: RHS
% sol_s: Solution
% Loop sum_level
%------------------------------------------------
dim=2;

% ['Generate 1D keys']
% tic
nx=[];px=[];kx=[];
for Lx=0:n
    for Px=0:2^max(0,Lx-1)-1
        for Kx=1:k
            nx=[nx;Lx];
            px=[px;Px];
            kx=[kx;Kx];
            
        end
    end
end
%
Key1dMesh=[nx,px,kx];
% toc

['Generate ',num2str(dim),'-D Hash Table']
% tic
Hash =Mesh_interface;
Hash.init (dim,n,k-1,k^dim);
% toc
dof_sparse = double(Hash.size);

% sort the index
All_index=allv(Hash);
All_index=sort(All_index);
Real_index=All_index'+1;

A_s = sparse(dof_sparse,dof_sparse);


b_s = sparse(dof_sparse,1);
sol_s = sparse(dof_sparse,1);
uu_s=sparse(dof_sparse,1);

count=1;
% Method 2
for sum_level=0:n
    for ix_level=0:sum_level
        iy_level=sum_level-ix_level;
        
        Ix=Index_1D(k,ix_level);
        Iy=Index_1D(k,iy_level);
        
        key_i=GenerateKey(Ix,Iy,Key1dMesh);
        Index_I=findvs(Hash,key_i);
        
        % Term 1: S*I+I*S at the diagonal entry
        tmp=kron(Stiff_1D(Ix,Ix),M_mass(Iy,Iy))+kron(M_mass(Ix,Ix),Stiff_1D(Iy,Iy));
        [II,JJ]=meshgrid(Real_index(Index_I+1),Real_index(Index_I+1));

        A_s=A_s+sparse(double(II),double(JJ),tmp,dof_sparse,dof_sparse);
        
        tmp=kron(b(Ix),b(Iy));

        b_s(Real_index(Index_I+1))=b_s(Real_index(Index_I+1))+tmp;
        uu_s(Real_index(Index_I+1))=uu_s(Real_index(Index_I+1))+kron(coef_MW(Ix),coef_MW(Iy));
        
                % save matrices to files
        A_encode{count}.A1=Stiff_1D(Ix,Ix);
        A_encode{count}.A2=M_mass(Iy,Iy);
        A_encode{count}.B1=M_mass(Ix,Ix);
        A_encode{count}.B2=Stiff_1D(Iy,Iy);
        A_encode{count}.IndexI=double(Real_index(Index_I+1));
        A_encode{count}.IndexJ=double(Real_index(Index_I+1));
        count=count+1;
        
        % Term 2: SxI--Assume jy_level==iy_level
        jy_level=iy_level;
        for jx_level=0:ix_level-1
            
            Jx=Index_1D(k,jx_level);
            Jy=Iy;

            key_j=GenerateKey(Jx,Jy,Key1dMesh);
            
            Index_J=findvs(Hash,key_j);
            
            [II,JJ]=meshgrid(Real_index(Index_J+1),Real_index(Index_I+1));
            
            
            
            tmp=kron(Stiff_1D(Jx,Ix),M_mass(Jy,Iy));
            A_s=A_s+sparse([double(II),double(JJ)],[double(JJ),double(II)],[tmp',tmp'],dof_sparse,dof_sparse);
            
            % save matrices to files
            A_encode{count}.A1=Stiff_1D(Jx,Ix);
            A_encode{count}.A2=M_mass(Jy,Iy);
            A_encode{count}.B1=0;
            A_encode{count}.B2=0;
            A_encode{count}.IndexI=double(Real_index(Index_J+1));
            A_encode{count}.IndexJ=double(Real_index(Index_I+1));

            A_encode{count+1}.A1=Stiff_1D(Jx,Ix)';
            A_encode{count+1}.A2=M_mass(Jy,Iy)';
            A_encode{count+1}.B1=0;
            A_encode{count+1}.B2=0;
            A_encode{count+1}.IndexI=double(Real_index(Index_I+1));
            A_encode{count+1}.IndexJ=double(Real_index(Index_J+1));
            
            count=count+2;
            
        end
        % Term 3: S*I
        jx_level=ix_level;
        for jy_level=0:iy_level-1
            Jx=Ix;
            Jy=Index_1D(k,jy_level);
     
            key_j=GenerateKey(Jx,Jy,Key1dMesh);
            
            Index_J=findvs(Hash,key_j);
            
            [II,JJ]=meshgrid(Real_index(Index_J+1),Real_index(Index_I+1));
            
            
            tmp=kron(M_mass(Jx,Ix),Stiff_1D(Jy,Iy));
            
            A_s=A_s+sparse([double(II),double(JJ)],[double(JJ),double(II)],[tmp',tmp'],dof_sparse,dof_sparse);
            
                        % save matrices to files
            A_encode{count}.A1=0;
            A_encode{count}.A2=0;
            A_encode{count}.B1=M_mass(Jx,Ix);
            A_encode{count}.B2=Stiff_1D(Jy,Iy);
            A_encode{count}.IndexI=double(Real_index(Index_J+1));
            A_encode{count}.IndexJ=double(Real_index(Index_I+1));

            A_encode{count+1}.A=0;
            A_encode{count+1}.B=0;
            A_encode{count+1}.C=M_mass(Jx,Ix)';
            A_encode{count+1}.D=Stiff_1D(Jy,Iy)';
            A_encode{count+1}.IndexI=double(Real_index(Index_I+1));
            A_encode{count+1}.IndexJ=double(Real_index(Index_J+1));
            
            count=count+2;
            
        end
        
          
    end
    
end
figure;
spy(A_s)
title(sprintf('4D problem, n=%d,nnz=%g,condest=%g',...
    size(A_s,1),nnz(A_s),condest(A_s)  ));
% Check matrix
% eigs(A_s,3,'SM')

save(['./Data/A_2D_encode.mat'],'A_encode');

% tic
sol_s = A_s\b_s*pi^2*2;
% toc

% compare the solution with interpolation
norm(sol_s-uu_s)
function Ix=Index_1D(k,level)

if level==0
    Ix=[1:k];
else
    Ix=[k*2^(level-1)+1:k*2^level];
end

end

function key=GenerateKey(Ix,Iy,Key1dMesh)

tmp_x=Key1dMesh(Ix,:);
tmp_y=Key1dMesh(Iy,:);
tmp_ix=reshape(repmat(tmp_x,1,size(tmp_y,1))',3,size(tmp_x,1)*size(tmp_y,1) )';
tmp_iy=repmat(tmp_y,size(tmp_x,1),1);
key=[tmp_ix(:,1),tmp_iy(:,1),tmp_ix(:,2),tmp_iy(:,2),tmp_ix(:,3)-1,tmp_iy(:,3)-1];

end
