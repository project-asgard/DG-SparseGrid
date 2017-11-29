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
run hashTest2D.m

A_s = sparse(dof_sparse,dof_sparse);
b_s = sparse(dof_sparse,1);
sol_s = sparse(dof_sparse,1);
uu_s=sparse(dof_sparse,1);

kron_flops = 0;
kron_nnz = 0;

count=1;
% Method 2
for sum_level=0:n
    for ix_level=0:sum_level
        iy_level=sum_level-ix_level;
        
        Ix=Index_1D(k,ix_level);
        Iy=Index_1D(k,iy_level);
        
        key_i=GenerateKey2D(Ix,Iy,Key1dMesh);
        for iii=1:size(key_i,1)
            Index_I(iii,1)=database.(sprintf('i%g_',key_i(iii,:)));
        end
        
        % Term 1: S*I+I*S at the diagonal entry
        tmp=kron(Stiff_1D(Ix,Ix),M_mass(Iy,Iy))+kron(M_mass(Ix,Ix),Stiff_1D(Iy,Iy));
        [II,JJ]=meshgrid(Index_I,Index_I);
        
        A_s=A_s+sparse(II,JJ,tmp,dof_sparse,dof_sparse);
        
        tmp=kron(b(Ix),b(Iy));
        
        b_s(Index_I)=b_s(Index_I)+tmp;
        uu_s(Index_I)=uu_s(Index_I)+kron(coef_MW(Ix),coef_MW(Iy));
        
        % save matrices to files
        A_encode{count}.A1=Stiff_1D(Ix,Ix);
        A_encode{count}.A2=M_mass(Iy,Iy);
        A_encode{count}.B1=M_mass(Ix,Ix);
        A_encode{count}.B2=Stiff_1D(Iy,Iy);
        A_encode{count}.IndexI=double(Index_I);
        A_encode{count}.IndexJ=double(Index_I);

        kron_flops = kron_flops + ...
                     kron_mult_cost2( ...
                        A_encode{count}.A1, ...
                        A_encode{count}.A2);

        kron_flops = kron_flops + ...
                     kron_mult_cost2( ...
                        A_encode{count}.B1, ...
                        A_encode{count}.B2);

        kron_nnz = kron_nnz + ...
                     length(A_encode{count}.IndexI ) + ...
                     length(A_encode{count}.IndexJ );
                      
                        
                     
               

        count=count+1;
        
        % Term 2: SxI--Assume jy_level==iy_level
        jy_level=iy_level;
        for jx_level=0:ix_level-1
            
            Jx=Index_1D(k,jx_level);
            Jy=Iy;
            
            key_j=GenerateKey2D(Jx,Jy,Key1dMesh);
            for jjj=1:size(key_j,1)
                Index_J(jjj,1)=database.(sprintf('i%g_',key_j(jjj,:)));
            end
            
            [II,JJ]=meshgrid(Index_J,Index_I);
            
            tmp=kron(Stiff_1D(Jx,Ix),M_mass(Jy,Iy));
            A_s=A_s+sparse([II,JJ],[JJ,II],[tmp',tmp'],dof_sparse,dof_sparse);
            
            % save matrices to files
            A_encode{count}.A1=Stiff_1D(Jx,Ix);
            A_encode{count}.A2=M_mass(Jy,Iy);
            A_encode{count}.B1=0;
            A_encode{count}.B2=0;
            A_encode{count}.IndexI=Index_J;
            A_encode{count}.IndexJ=Index_I;

            kron_flops = kron_flops + ...
                     kron_mult_cost2( ...
                        A_encode{count}.A1, ...
                        A_encode{count}.A2);


            kron_nnz = kron_nnz + ...
                     length(A_encode{count}.IndexI ) + ...
                     length(A_encode{count}.IndexJ );
                      
            
            A_encode{count+1}.A1=Stiff_1D(Jx,Ix)';
            A_encode{count+1}.A2=M_mass(Jy,Iy)';
            A_encode{count+1}.B1=0;
            A_encode{count+1}.B2=0;
            A_encode{count+1}.IndexI=Index_I;
            A_encode{count+1}.IndexJ=Index_J;
            
            kron_flops = kron_flops + ...
                     kron_mult_cost2( ...
                        A_encode{count+1}.A1, ...
                        A_encode{count+1}.A2);


            kron_nnz = kron_nnz + ...
                     length(A_encode{count+1}.IndexI ) + ...
                     length(A_encode{count+1}.IndexJ );



            count=count+2;
            clear Index_J
        end
        % Term 3: S*I
        jx_level=ix_level;
        for jy_level=0:iy_level-1
            Jx=Ix;
            Jy=Index_1D(k,jy_level);
            
            key_j=GenerateKey2D(Jx,Jy,Key1dMesh);
            for jjj=1:size(key_j,1)
                Index_J(jjj,1)=database.(sprintf('i%g_',key_j(jjj,:)));
            end
            
            
            [II,JJ]=meshgrid(Index_J,Index_I);
            
            
            tmp=kron(M_mass(Jx,Ix),Stiff_1D(Jy,Iy));
            
            A_s=A_s+sparse([II,JJ],[JJ,II],[tmp',tmp'],dof_sparse,dof_sparse);
            
            % save matrices to files
            A_encode{count}.A1=0;
            A_encode{count}.A2=0;
            A_encode{count}.B1=M_mass(Jx,Ix);
            A_encode{count}.B2=Stiff_1D(Jy,Iy);
            A_encode{count}.IndexI=Index_J;
            A_encode{count}.IndexJ=Index_I;

            kron_flops = kron_flops + ...
                     kron_mult_cost2( ...
                        A_encode{count}.B1, ...
                        A_encode{count}.B2);


            kron_nnz = kron_nnz + ...
                     length(A_encode{count}.IndexI ) + ...
                     length(A_encode{count}.IndexJ );
            
%EFD  ---------------------------------------
%EFD  use A1,A2,B1,B2 instead of A,B,C,D
%EFD 
%EFD             A_encode{count+1}.A=0;
%EFD             A_encode{count+1}.B=0;
%EFD             A_encode{count+1}.C=M_mass(Jx,Ix)';
%EFD             A_encode{count+1}.D=Stiff_1D(Jy,Iy)';
%EFD  ------------------------------------

            A_encode{count+1}.A1=0;
            A_encode{count+1}.A2=0;
            A_encode{count+1}.B1=M_mass(Jx,Ix)';
            A_encode{count+1}.B2=Stiff_1D(Jy,Iy)';
            A_encode{count+1}.IndexI=Index_I;
            A_encode{count+1}.IndexJ=Index_J;

            kron_flops = kron_flops + ...
                     kron_mult_cost2( ...
                        A_encode{count+1}.B1, ...
                        A_encode{count+1}.B2);


            kron_nnz = kron_nnz + ...
                     length(A_encode{count+1}.IndexI ) + ...
                     length(A_encode{count+1}.IndexJ );
            
            count=count+2;
            clear Index_J
        end
        
        clear Index_I
    end
    
end
% -----------------------
% ensure A_s is symmetric
% -----------------------
A_s = (A_s + A_s')/2;

disp(sprintf('Np=%d,k=%d, kron_flops=%g, kron_nnz=%g, nnz(A_s)=%g', ...
              Np, k, kron_flops,    kron_nnz, nnz(A_s)  ));

figure;
spy(A_s)
title(sprintf('2D problem, Np=%d,k=%d, n=%d,nnz=%g,condest=%g',...
    Np, k, ...
    size(A_s,1),nnz(A_s),condest(A_s)  ));
% Check matrix
% -------------------------
% estimate condition number
% -------------------------
% eigs( A_s, 2, 'BE');

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

function key=GenerateKey2D(Ix,Iy,Key1dMesh)

tmp_x=Key1dMesh(Ix,:);
tmp_y=Key1dMesh(Iy,:);
tmp_ix=reshape(repmat(tmp_x,1,size(tmp_y,1))',3,size(tmp_x,1)*size(tmp_y,1) )';
tmp_iy=repmat(tmp_y,size(tmp_x,1),1);
key=[tmp_ix(:,1),tmp_iy(:,1),tmp_ix(:,2),tmp_iy(:,2),tmp_ix(:,3)-1,tmp_iy(:,3)-1];

end
