function [flops1, flops2, imethod] = flops_kron3( nrow1,ncol1, ...
                                                  nrow2,ncol2, ...
                                                  nrow3,ncol3)

% [flops1, flops2, imethod] = flops_kron3( nrow1,ncol1, ...
%                                          nrow2,ncol2, ...
%                                          nrow3,ncol3)
%
% estimate flops for method 1 and method 2
%

% ----------------------
% method 1:
% Ytmp = kron(A2,A3)*X
% Y = Ytmp * transpose(A1)
% ----------------------
nrowX = ncol1*ncol2*ncol3;
nvec = ncol1;
flops1 = flops_kron2( nrow2,ncol2, nrow3,ncol3)*nvec;
%
% Ytmp is (nrow2*nrow3) by ncol1;
% 
nrow_Ytmp = (nrow2*nrow3);  
flops1 = flops1 + 2.0 * nrow_Ytmp * ncol1 * nrow1;


% -------------- 
% method 2:
% Ytmp =   X * transpose(A1);
% Y = kron(A2,A3)*Ytmp
% -------------- 

% -------------------------------------
% X reshaped as  (ncol2*ncol3) by ncol1
% Ytmp is (ncol2*ncol3) by nrow1
% -------------------------------------
flops2 =  2.0 * (ncol2*ncol3) * ncol1 * nrow1;

flops2 = flops2 + flops_kron2( nrow2,ncol2, nrow3,ncol3)*nrow1;

if (flops1 <= flops2),
  imethod = 1;
else
  imethod = 2;
end;
