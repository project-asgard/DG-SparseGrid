function [fnew]=multi_2D(A,B,f,HASHInv,Lev,Deg)
%====================================================================
% This code is to compute the product of Kron(A,B)*f
% Input: Matrix: A (size nA1*nA2)
%        Matrix: B (size nB1*nB2)
%        Vector: f (size nf*1)
%        HASHInv: inverse of the Hash Table
% Output: Vector: fnew (size nf*1)
%====================================================================


% Sparse Grid DG Basis
HASHDOF = size(HASHInv,2);

dof_1D_FG=Deg*2^(Lev);

% output a new vector with full-grid size
fnew=sparse(dof_1D_FG^2,1);

for ii=1:HASHDOF
    
    I1=HASHInv{ii}(5);
%     n1=HashInv{ii}(1);
    I2=HASHInv{ii}(6);
%     n2=HashInv{ii}(2);
    
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



