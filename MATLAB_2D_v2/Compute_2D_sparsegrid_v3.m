% Computing on the sparse grid for Poisson Eq
% Key1DMesh
['Generate 1D keys']
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

['Generate ',num2str(2),'-D Hash Table']
% tic
Hash =Mesh_interface;
Hash.init (2,n,k-1,20);
% toc
dof_sparse = double(Hash.size);

% sort the index
All_index=allv(Hash);
All_index=sort(All_index);
Real_index=All_index'+1;

A_s = sparse(dof_sparse,dof_sparse);


b_s = sparse(dof_sparse,1);
sol_s = sparse(dof_sparse,1);


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
        b_s(Index_I+1)=b_s(Index_I+1)+tmp;
              
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
            
            
        end
        
          
    end
    
end

spy(A_s)

eigs(A_s,3,'SM')

tic
sol_s = A_s\b_s*pi^2*2;
toc
% full(sol_s)

% % Method 3 not work
% dim=2;n=2;
% for sum_level=0:n
% %     sum_level
%     n_level=DFS(dim,sum_level,0,[],[]);
% %     level
%     
%     for i=1:size(n_level,1)
%         n_level(i,:)
%         % Diagonal Entry 
%         
%         
%         for dd=1:dim
%         Id=Index_1D(k,n_level(i,dd));
%         Jd=Index_1D(k,m_level(i,dd));
%         
%         tmp1=kron(Stiff_1D(Id,Id),M_mass(Jy,Iy));
%         tmp2=kron(Stiff_1D(Jx,Ix),M_mass(Jy,Iy));
%         end
%         
% %         key_i=GenerateKey(Id,Key1dMesh);
%         
% %         I1=Index_1D(k,n_level(i,1));
% %         I2=Index_1D(k,n_level(i,2));
% %         
% %         key_i=GenerateKey(I1,I2,Key1dMesh);
% %         Index_I=findvs(Hash,key_i);
%         
%         m_level=nestedLoopOperation(zeros(1,dim), n_level(i,:), 1);
%         
%         
%     end
% end

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

% function key=GenerateKeyd(Id,Key1dMesh)
% % for d-dimension
% for dd=1:size(Id,1)
% tmp_d=Key1dMesh(Id(dd),:);
% end
% 
% tmp_y=Key1dMesh(Iy,:);
% tmp_ix=reshape(repmat(tmp_x,1,size(tmp_y,1))',3,size(tmp_x,1)*size(tmp_y,1) )';
% tmp_iy=repmat(tmp_y,size(tmp_x,1),1);
% key=[tmp_ix(:,1),tmp_iy(:,1),tmp_ix(:,2),tmp_iy(:,2),tmp_ix(:,3)-1,tmp_iy(:,3)-1];
% 
% end


% return
% method 2: Working!
% count=1;
% % Stiffness Matrix
% for my=0:n
%     if my==0
%         Iy=[1:k];
%         IndexI=[1:k^2*2^n];
%     else
%         Iy=[k*2^(my-1)+1:k*2^my];
%         IndexI=[(my-1)*k^2*2^(n-1)+k^2*2^n+1:(my)*k^2*2^(n-1)+k^2*2^n];
%     end
%
%     mx=[0:n-my];
%     Ix=[1:k*2^(n-my)];
%
%     for ny=0:n
%
%         if ny==0
%             Jy=[1:k];
%             IndexJ=[1:k^2*2^n];
%         else
%             Jy=[k*2^(ny-1)+1:k*2^ny];
%             IndexJ=[(ny-1)*k^2*2^(n-1)+k^2*2^n+1:(ny)*k^2*2^(n-1)+k^2*2^n];
%         end
%
%         nx=[0:n-ny];
%         Jx=[1:k*2^(n-ny)];
%
%         tmp=(kron(M_mass(Iy,Jy),Stiff_1D(Ix,Jx))+kron(Stiff_1D(Iy,Jy),M_mass(Ix,Jx)));
%
%         % save matrices to files
%         Mass_tmp=M_mass(Iy,Jy);Stiff_tmp=Stiff_1D(Ix,Jx);
%         save(['./Data/M_mass1_',num2str(count),'.mat'],'Mass_tmp','Stiff_tmp')
%         Mass_tmp=M_mass(Ix,Jx);Stiff_tmp=Stiff_1D(Iy,Jy);
%         save(['./Data/M_mass2_',num2str(count),'.mat'],'Mass_tmp','Stiff_tmp')
%         save(['./Data/Index_',num2str(count),'.mat'],'IndexI','IndexJ')
%         count=count+1;
%
%         [xindex,yindex]=meshgrid(IndexI,IndexJ);
%
%         A_s=A_s+sparse(xindex',yindex',tmp,dof_sparse,dof_sparse);
%
%
%     end
%     tmp=kron(b(Iy),b(Ix));
%     b_s(IndexI)=b_s(IndexI)+tmp;
%
% end
%
% tic
% sol_s = A_s\b_s*pi^2*2;
% toc
%
% cond_sparse = condest(A_s);
%
% % generate S2F matrix
% S2F = sparse(dof_sparse,dof_full);
% for my=0:n
%     % position for sparse grid
%     %     my
%     if my==0
%
%         IndexI=[1:k^2*2^n];
%         Iv=[1:k^2*2^n];
%     else
%
%         IndexI=[(my-1)*k^2*2^(n-1)+k^2*2^n+1:(my)*k^2*2^(n-1)+k^2*2^n];
%         Iv=[(k*2^(my-1))*(k*2^n)+1:(k*2^(my))*(k*2^n)];
%     end
%
%
%     tmp=[];
%     % delete entry from Iv
%     for iy=0:2^max(0,my-1)-1
%         mx=n-my;
%         for indexk=1:k
%             if mx<n
%                 tmp=[tmp,...
%                     iy*(k^2*2^(n))+(indexk-1)*(k*2^n)+[(k*2^(mx)+1:k*2^(n))],...
%                     ];
%             end
%         end
%     end
%
%     Iv(tmp)=[];
%
%     S2F=S2F+sparse(IndexI,Iv,ones(1,size(IndexI,2)),dof_sparse,dof_full);
%
% end
%
%
% Lmax_2D_s = max(abs(Iu_2D-S2F'*sol_s));