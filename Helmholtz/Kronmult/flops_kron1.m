function [flops1] = flops_kron1( nrow1, ncol1, nvec )
% [flops1] = flops_kron1( nrow1, ncol1, nvec )
%
%  flops for  A1 * X
%
flops1 = 2.0*nrow1*ncol1*nvec;

