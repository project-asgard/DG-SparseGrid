function [fnew]=Multi_2D(A,B,f,Hash,HashInv,Deg)
%====================================================================
% This code is to compute the product of Kron(A,B)*f
% Input: Matrix: A (size nA1*nA2)
%        Matrix: B (size nB1*nB2)
%        Vector: f (size nf*1)
%        Hash Table: Hash
%        isFull: isFull=1--output the fullgrid vector
%                    isFull=0--output the sparsegrid vector
% Output: Vector: fnew (size nf*1)
%====================================================================


% Sparse Grid DG Basis
Lev_v = Hash.Lev;
Lev_x = Hash.Lev;
Lev = Hash.Lev;

% Polynomial Degree


dof_1D_x=Deg*2^(Lev_x);
dof_1D_v=Deg*2^(Lev_v);

% output a new vector with full-grid size
fnew=sparse(dof_1D_v*dof_1D_x,1);

tic

for ii=1:Hash.dof
    %         ii
    %         n1=ll(1);c1=ll(3);
    %         n2=ll(2);c2=ll(4);
    
    I1=HashInv.x1(ii);%LevCell2index(n1,c1);
    n1=ceil(log2(I1));
    I2=HashInv.x2(ii);%LevCell2index(n2,c2);
    n2=ceil(log2(I2));
    
    index_I1=[(I1-1)*Deg+1:I1*Deg];
    index_I2=[(I2-1)*Deg+1:I2*Deg];
    
    Index = Deg^2*(ii-1)+1:Deg^2*ii;
    %         [n1 n2 Index]
    
    tmp=kron(...
        A(:,index_I1),...
        B(:,index_I2) ...
        )*f(Index(:));
    fnew=fnew+tmp;
    
    
    %         Deg*(HashInv.x1(ii)-1)+1:Deg*HashInv.x1(ii)
    %         Deg*(HashInv.x2(ii)-1)+1:Deg*HashInv.x2(ii)
    %         2^n1*2^n2
    % %         +[Deg*(ii-1)+1:Deg*ii]
    % %         Deg*(HashInv.x2(ii)-1)+1:Deg*HashInv.x2(ii)
    %         disp('========')
    % % size(A(:,Deg*(HashInv.x1(ii)-1)+1:Deg*HashInv.x1(ii)))
    % % size(B(:,Deg*(HashInv.x2(ii)-1)+1:Deg*HashInv.x2(ii)))
    % % Deg*(HashInv.x1(ii)-1)+1:Deg*HashInv.x1(ii))
    %
    % % size(kron(...
    % %                             A(:,Deg*(HashInv.x1(ii)-1)+1:Deg*HashInv.x1(ii)),...
    % %                             B(:,Deg*(HashInv.x2(ii)-1)+1:Deg*HashInv.x2(ii)) ...
    % %                            ))
    % %   size(f(Deg^2*(ii-1)+1:Deg^2*ii))
    %         tmp=kron(...
    %                             A(:,Deg*(HashInv.x1(ii)-1)+1:Deg*HashInv.x1(ii)),...
    %                             B(:,Deg*(HashInv.x2(ii)-1)+1:Deg*HashInv.x2(ii)) ...
    %                            )*f(Deg^2*(ii-1)+1:Deg^2*ii);
    %         fnew=fnew+tmp;
end


end



