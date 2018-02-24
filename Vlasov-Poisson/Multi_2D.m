function [fnew]=Multi_2D(A,B,f,Hash,isFull,HashInv)
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
if nargin == 4
    isFull=1;
end

% Sparse Grid DG Basis
Lev_v = Hash.Lev_1;
Lev_x = Hash.Lev_2;

% Polynomial Degree
Deg = Hash.Deg;

dof_1D_x=Deg*2^(Lev_x);
dof_1D_v=Deg*2^(Lev_v);
if isFull ==1
    % output a new vector with full-grid size
    fnew=sparse(dof_1D_v*dof_1D_x,1);
    
    tic
    
    for ii=1:Hash.dof
        tmp=kron(A(:,HashInv.x1(ii)),B(:,HashInv.x2(ii)))*f(ii);
        fnew=fnew+tmp;
    end
    
elseif isFull == 0
    % output a vector with sparse grids size
    fnew=sparse(Hash.dof,1);
    tic
    for ii=1:Hash.dof
        I_v=HashInv.x1(ii);
        I_x=HashInv.x2(ii);
        tmp=A(I_v,HashInv.x1).*B(I_x,HashInv.x2);
        fnew(ii)=fnew(ii)+tmp*(f);
    end
    
    
end





end



