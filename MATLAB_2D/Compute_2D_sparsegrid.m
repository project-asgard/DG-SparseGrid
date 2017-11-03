% Computing on the sparse grid for Poisson Eq

dof_sparse = n*k^2*2^(n-1)+k^2*2^n;

A_s = sparse(dof_sparse,dof_sparse);
b_s = sparse(dof_sparse,1);
sol_s = sparse(dof_sparse,1);

count=1;
% Stiffness Matrix
for my=0:n
    if my==0
        Iy=[1:k];
        IndexI=[1:k^2*2^n];
    else
        Iy=[k*2^(my-1)+1:k*2^my];
        IndexI=[(my-1)*k^2*2^(n-1)+k^2*2^n+1:(my)*k^2*2^(n-1)+k^2*2^n];
    end
    
    mx=[0:n-my];
    Ix=[1:k*2^(n-my)];
    
    for ny=0:n
        
        if ny==0
            Jy=[1:k];
            IndexJ=[1:k^2*2^n];
        else
            Jy=[k*2^(ny-1)+1:k*2^ny];
            IndexJ=[(ny-1)*k^2*2^(n-1)+k^2*2^n+1:(ny)*k^2*2^(n-1)+k^2*2^n];
        end
        
        nx=[0:n-ny];
        Jx=[1:k*2^(n-ny)];
        
        tmp=(kron(M_mass(Iy,Jy),Stiff_1D(Ix,Jx))...
            +kron(Stiff_1D(Iy,Jy),M_mass(Ix,Jx)));
        
        % save matrices to files
        A_encode{count}.A=M_mass(Iy,Jy);
        A_encode{count}.B=Stiff_1D(Ix,Jx);
        A_encode{count}.C=Stiff_1D(Iy,Jy);
        A_encode{count}.D=M_mass(Ix,Jx);
        A_encode{count}.IndexI=IndexI;
        A_encode{count}.IndexJ=IndexJ;
        
        
        count=count+1;
        
        [xindex,yindex]=meshgrid(IndexI,IndexJ);
        
        A_s=A_s+sparse(xindex',yindex',tmp,dof_sparse,dof_sparse);
        
        
    end
    tmp=kron(b(Iy),b(Ix));
    b_s(IndexI)=b_s(IndexI)+tmp;
    
end

save(['./Data/A_encode.mat'],'A_encode');

tic
sol_s = A_s\b_s*pi^2*2;
toc

cond_sparse = condest(A_s);

% generate S2F matrix 
S2F = sparse(dof_sparse,dof_full);
for my=0:n
    % position for sparse grid
    %     my
    if my==0
        
        IndexI=[1:k^2*2^n];
        Iv=[1:k^2*2^n];
    else
        
        IndexI=[(my-1)*k^2*2^(n-1)+k^2*2^n+1:(my)*k^2*2^(n-1)+k^2*2^n];
        Iv=[(k*2^(my-1))*(k*2^n)+1:(k*2^(my))*(k*2^n)];
    end
    
    
    tmp=[];
    % delete entry from Iv
    for iy=0:2^max(0,my-1)-1
        mx=n-my;
        for indexk=1:k
            if mx<n
                tmp=[tmp,...
                    iy*(k^2*2^(n))+(indexk-1)*(k*2^n)+[(k*2^(mx)+1:k*2^(n))],...
                    ];
            end
        end
    end
    
    Iv(tmp)=[];
    
    S2F=S2F+sparse(IndexI,Iv,ones(1,size(IndexI,2)),dof_sparse,dof_full);
    
end


Lmax_2D_s = max(abs(Iu_2D-S2F'*sol_s));