function [fnew]=Multi_2D(A,B,f,Hash,HashInv,Deg)
%====================================================================
% This code is to compute the product of Kron(A,B)*f
% Input: Matrix: A (size nA1*nA2)
%        Matrix: B (size nB1*nB2)
%        Vector: f (size nf*1)
%        Hash Table: Hash
%        HashInv: inverse of the Hash Table
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
 
    I1=HashInv{ii}(5);
    n1=HashInv{ii}(1);
    I2=HashInv{ii}(6);
    n2=HashInv{ii}(2);
    
    index_I1=[(I1-1)*Deg+1:I1*Deg];
    index_I2=[(I2-1)*Deg+1:I2*Deg];
    
    Index = Deg^2*(ii-1)+1:Deg^2*ii;
    
    tmp=kron(...
        A(:,index_I1),...
        B(:,index_I2) ...
        )*f(Index(:));
    fnew=fnew+tmp;

end


end



