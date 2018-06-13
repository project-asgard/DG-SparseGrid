function [flops1, flops2, imethod] = flops_kron2( nrow1,ncol1, ...
                                                  nrow2,ncol2, ...
                                                  nvec)

% [flops1, flops2, imethod] = flops_kron2( nrow1,ncol1, ...
%                                                   nrow2,ncol2, ...
%                                                   nvec)
%
% estimate flops for method 1 and method 2
%

% ----------------------
% method 1:
% A2X = A2*X
% Y = A2X * transpose(A1)
% ----------------------
flops1 = 2.0*nrow2*ncol2*nvec;
nrow_A2X = nrow2;
flops1 = flops1 + 2.0 * nrow_A2X * ncol1 * nrow1 ;


% -------------- 
% method 2:
% XA1t = X * transpose(A1)
% Y = A2 * XA1t
% -------------- 
nrow_X = ncol2;
flops2 = 2.0 * nrow_X * ncol1 * nrow1;

ncol_XA1t = nrow1;
flops2 = flops2 + 2.0 * nrow2 * ncol2 * ncol_XAt1;

if (flops1 <= flops2),
  imethod = 1;
else
  imethod = 2;
end;
