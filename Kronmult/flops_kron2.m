function [flops1, flops2, imethod, imem1, imem2] = flops_kron2( nrow1,ncol1, ...
                                                  nrow2,ncol2)

% [flops1, flops2, imethod, imem1, imem2] = flops_kron2( nrow1,ncol1, ...
%                                          nrow2,ncol2 )
%
% estimate flops (and temporary storage) for method 1 and method 2
%

% ----------------------
% method 1:
% A2X = A2*X
% Y = A2X * transpose(A1)
% ----------------------

% -------------------
% X is ncol2 by ncol1
% -------------------
flops1 = 2.0*nrow2*ncol2 * ncol1;
nrow_A2X = nrow2;
ncol_A2X = nrow1;
imem1 = nrow_A2X * ncol_A2X;
flops1 = flops1 + 2.0 * nrow_A2X * ncol1 * nrow1 ;


% -------------- 
% method 2:
% XA1t = X * transpose(A1)
% Y = A2 * XA1t
% -------------- 
nrow_X = ncol2;
flops2 = 2.0 * nrow_X * ncol1 * nrow1;

ncol_XA1t = nrow1;
nrow_XA1t = nrow_X;
imem2 = nrow_XA1t * ncol_XA1t;
flops2 = flops2 + 2.0 * nrow2 * ncol2 * ncol_XA1t;

if (flops1 <= flops2),
  imethod = 1;
else
  imethod = 2;
end;
