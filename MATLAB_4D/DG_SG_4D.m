% Test 4D problem
% (n1,i1,p1)
% (n2,i2,p2)
% (n3,i3,p3)
% (n4,i4,p4)
% Just check the stiffness matrix for Poisson equation
%--------------------------------------------------------------
clc
clear all
close all
format short e

%% 1. Set up parameters
n=5;k=2;
Dim=4;

%% 2. Generate Hash Table for 1D
['Generate 1D Hash Table']
tic
% initial meshing
Mesh1D = Mesh_interface;
Mesh1D.init(1,n,k-1,12);
toc

%% 3. Generate Hash table for 4D
['Generate ',num2str(Dim),'-D Hash Table']
tic
Hash =Mesh_interface;
Hash.init (Dim,n,k-1,20);
toc
dof_sparse = double(Hash.size);

% sort the index
All_index=allv(Hash);
All_index=sort(All_index);
Real_index=All_index'+1;

%% 4. 1D matrix generation
['Generate 1D Laplacian Matrix']
tic
[Stiff_1D,b] = LaplacianMatrix2(n,k);
% [Stiff_1D,b1,b2,xnode,Meval,uu]=LaplacianMatrix(N_max,k);
% load(['Gauss_L12_k',num2str(k),'.mat'])
toc

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

%% 5. Construct 4D stiffness matrix
A_s = sparse(dof_sparse,dof_sparse);
b_s = sparse(dof_sparse,1);
sol_s = sparse(dof_sparse,1);

['Construct Matrix']
tic % varying n1
for n2=0:n
    for p2=0:2^max(n2-1,0)-1
        for k2=1:k % for Group1
            
            for n3=0:max(n-n2,0)
                for p3=0:2^max(n3-1,0)-1
                    for k3=1:k % for Group2
                        
                        for n4=0:max(n-n2-n3,0)
                            for p4=0:2^max(n4-1,0)-1
                                for k4=1:k % for Group3
                                    
                                    tmp=Key1dMesh(1:2^((n-n2-n3-n4))*k,:);
                                    
                                    n1=tmp(:,1);
                                    p1=tmp(:,2);
                                    k1=tmp(:,3);
                                    
                                    zz=ones(size(n1,1),1);
                                    keys=[n1,n2*zz,n3*zz,n4*zz,...
                                        p1,p2*zz,p3*zz,p4*zz,...
                                        (k1-1),(k2-1)*zz,(k3-1)*zz,(k4-1)*zz];
                                    val=Stiff_1D(1:size(n1,1),1:size(n1,1));
                                    Index_J=findvs(Hash,keys);
                                    
                                    [II,JJ]=meshgrid(Real_index(Index_J+1));
                                    A_s=A_s+sparse(double(II),double(JJ),val,dof_sparse,dof_sparse);
                                    
                                end
                            end
                        end
                        
                        
                    end
                end
            end
            
            
        end
    end
end
toc


tic 
for n1=0:n
    for p1=0:2^max(n1-1,0)-1
        for k1=1:k % for Group1
            
            for n3=0:max(n-n1,0)
                for p3=0:2^max(n3-1,0)-1
                    for k3=1:k % for Group2
                        
                        for n4=0:max(n-n1-n3,0)
                            for p4=0:2^max(n4-1,0)-1
                                for k4=1:k % for Group3
                                    
                                    tmp=Key1dMesh(1:2^((n-n1-n3-n4))*k,:);
                                    
                                    n2=tmp(:,1);
                                    p2=tmp(:,2);
                                    k2=tmp(:,3);
                                    
                                    zz=ones(size(n2,1),1);
                                    keys=[n1*zz,n2,n3*zz,n4*zz,...
                                        p1*zz,p2,p3*zz,p4*zz,...
                                        (k1-1)*zz,(k2-1),(k3-1)*zz,(k4-1)*zz];
                                    val=Stiff_1D(1:size(n2,1),1:size(n2,1)*zz);
                                    Index_J=findvs(Hash,keys);
                                    
                                    [II,JJ]=meshgrid(Real_index(Index_J+1));
                                    A_s=A_s+sparse(double(II),double(JJ),val,dof_sparse,dof_sparse);
                                    
                                end
                            end
                        end
                        
                        
                    end
                end
            end
            
            
        end
    end
end
toc

tic % varying n3
for n1=0:n
    for p1=0:2^max(n1-1,0)-1
        for k1=1:k % for Group1
            
            for n2=0:max(n-n1,0)
                for p2=0:2^max(n2-1,0)-1
                    for k2=1:k % for Group2
                        
                        for n4=0:max(n-n1-n2,0)
                            for p4=0:2^max(n4-1,0)-1
                                for k4=1:k % for Group3
                                    
                                    tmp=Key1dMesh(1:2^((n-n1-n2-n4))*k,:);
                                    
                                    n3=tmp(:,1);
                                    p3=tmp(:,2);
                                    k3=tmp(:,3);
                                    
                                    zz=ones(size(n3,1),1);
                                    keys=[n1*zz,n2*zz,n3,n4*zz,...
                                        p1*zz,p2*zz,p3,p4*zz,...
                                        (k1-1)*zz,(k2-1)*zz,(k3-1),(k4-1)*zz];
                                    val=Stiff_1D(1:size(n3,1),1:size(n3,1));
                                    Index_J=findvs(Hash,keys);
                                    
                                    [II,JJ]=meshgrid(Real_index(Index_J+1));
                                    A_s=A_s+sparse(double(II),double(JJ),val,dof_sparse,dof_sparse);
                                    
                                end
                            end
                        end
                        
                        
                    end
                end
            end
            
            
        end
    end
end
toc

tic % varying n4
for n1=0:n
    for p1=0:2^max(n1-1,0)-1
        for k1=1:k % for Group1
            
            for n2=0:max(n-n1,0)
                for p2=0:2^max(n2-1,0)-1
                    for k2=1:k % for Group2
                        
                        for n3=0:max(n-n1-n2,0)
                            for p3=0:2^max(n3-1,0)-1
                                for k3=1:k % for Group3
                                    
                                    tmp=Key1dMesh(1:2^((n-n1-n2-n3))*k,:);
                                    
                                    n4=tmp(:,1);
                                    p4=tmp(:,2);
                                    k4=tmp(:,3);
                                    
                                    zz=ones(size(n4,1),1);
                                    keys=[n1*zz,n2*zz,n3*zz,n4,...
                                        p1*zz,p2*zz,p3*zz,p4,...
                                        (k1-1)*zz,(k2-1)*zz,(k3-1)*zz,(k4-1)];
                                    val=Stiff_1D(1:size(n4,1),1:size(n4,1));
                                    Index_J=findvs(Hash,keys);
                                    
                                    [II,JJ]=meshgrid(Real_index(Index_J+1));
                                    A_s=A_s+sparse(double(II),double(JJ),val,dof_sparse,dof_sparse);
                                    
                                end
                            end
                        end
                        
                        
                    end
                end
            end
            
            
        end
    end
end
toc

figure;
spy(A_s)
['Check matrix']
tic
d=eigs((A_s),3,'SM')
toc
[d(1) 4*pi^2 d(1)-4*pi^2]

