% Computing on the sparse grid for Poisson Eq
% Key1DMesh
['Generate 1D keys']
tic
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
toc

['Generate ',num2str(2),'-D Hash Table']
tic
Hash =Mesh_interface;
Hash.init (2,n,k-1,20);
toc
dof_sparse = double(Hash.size);

% sort the index
All_index=allv(Hash);
All_index=sort(All_index);
Real_index=All_index'+1;

% dof_sparse = n*k^2*2^(n-1)+k^2*2^n;

A_s = sparse(dof_sparse,dof_sparse);
b_s = sparse(dof_sparse,1);
sol_s = sparse(dof_sparse,1);

% Method 2
for sum_level=0:n
    for ix_level=0:sum_level
        iy_level=sum_level-ix_level;
        % Term 1: SxI
        jy_level=iy_level;
        for jx_level=0:sum_level-jy_level
            Ix=Index_1D(k,ix_level);
            Iy=Index_1D(k,iy_level);
            Jx=Index_1D(k,jx_level);
            Jy=Iy;%Index_1D(k,jy_level);
            
            tmp_x=Key1dMesh(Ix,:);
            tmp_y=Key1dMesh(Iy,:);
            tmp_ix=reshape(repmat(tmp_x,1,size(tmp_y,1))',3,size(tmp_x,1)*size(tmp_y,1) )';
            tmp_iy=repmat(tmp_y,size(tmp_x,1),1);
            tmp_x=Key1dMesh(Jx,:);
            tmp_y=Key1dMesh(Jy,:);
            tmp_jx=reshape(repmat(tmp_x,1,size(tmp_y,1))',3,size(tmp_x,1)*size(tmp_y,1) )';
            tmp_jy=repmat(tmp_y,size(tmp_x,1),1);
            
            
            key_i=[tmp_ix(:,1),tmp_iy(:,1),tmp_ix(:,2),tmp_iy(:,2),tmp_ix(:,3)-1,tmp_iy(:,3)-1];
            key_j=[tmp_jx(:,1),tmp_jy(:,1),tmp_jx(:,2),tmp_jy(:,2),tmp_jx(:,3)-1,tmp_jy(:,3)-1];
            
            
             Index_I=findvs(Hash,key_i);
             Index_J=findvs(Hash,key_j);
                                    
             [II,JJ]=meshgrid(Real_index(Index_J+1),Real_index(Index_I+1));
             
            
%             tmp=kron(Stiff_1D(Jx,Ix),M_mass(Iy,Jy));
            tmp=kron(M_mass(Iy,Jy),Stiff_1D(Jx,Ix));

           
            if ix_level~=jx_level
                A_s=A_s+sparse([double(II),double(JJ)],[double(JJ),double(II)],[tmp',tmp'],dof_sparse,dof_sparse);
            else
                 A_s=A_s+sparse(double(II),double(JJ),tmp,dof_sparse,dof_sparse);
            end

            
        end
        % Term 2: IxS
        jx_level=ix_level;
        for jy_level=0:sum_level-jx_level
            Ix=Index_1D(k,ix_level);
            Iy=Index_1D(k,iy_level);
            Jx=Ix;%Index_1D(k,jx_level);
            Jy=Index_1D(k,jy_level);
            
            tmp_x=Key1dMesh(Ix,:);
            tmp_y=Key1dMesh(Iy,:);
            tmp_ix=reshape(repmat(tmp_x,1,size(tmp_y,1))',3,size(tmp_x,1)*size(tmp_y,1) )';
            tmp_iy=repmat(tmp_y,size(tmp_x,1),1);
            tmp_x=Key1dMesh(Jx,:);
            tmp_y=Key1dMesh(Jy,:);
            tmp_jx=reshape(repmat(tmp_x,1,size(tmp_y,1))',3,size(tmp_x,1)*size(tmp_y,1) )';
            tmp_jy=repmat(tmp_y,size(tmp_x,1),1);
            
            
            key_i=[tmp_ix(:,1),tmp_iy(:,1),tmp_ix(:,2),tmp_iy(:,2),tmp_ix(:,3)-1,tmp_iy(:,3)-1];
            key_j=[tmp_jx(:,1),tmp_jy(:,1),tmp_jx(:,2),tmp_jy(:,2),tmp_jx(:,3)-1,tmp_jy(:,3)-1];
            
            
             Index_I=findvs(Hash,key_i);
             Index_J=findvs(Hash,key_j);
                                    
             [II,JJ]=meshgrid(Real_index(Index_J+1),Real_index(Index_I+1));
            
%             tmp=kron(M_mass(Jx,Ix),Stiff_1D(Jy,Iy));
            tmp=kron(Stiff_1D(Jy,Iy),M_mass(Jx,Ix));
            
            if iy_level~=jy_level
                A_s=A_s+sparse([double(II),double(JJ)],[double(JJ),double(II)],[tmp',tmp'],dof_sparse,dof_sparse);
            else
                A_s=A_s+sparse(double(II),double(JJ),tmp,dof_sparse,dof_sparse);
            end
%             full(tmp)

            
        end
        
    end
    
end

spy(A_s)

eigs(A_s,3,'SM')
% Method 1


function Ix=Index_1D(k,level)
    if level==0
        Ix=[1:k];
    else
        Ix=[k*2^(level-1)+1:k*2^level];
    end

end


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