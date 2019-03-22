function Q = ortho( F )
% Q = ortho( F )
% generate ortho-normal basis from columns of matrix F
%
% essentially Gram-Schmidt (GS) orthogonalization but
% GS is numerically not as stable compared to QR factorization
% thus use QR instead
%

% ---------------------------------
% Note use "economy" version of QR
% so matrix Q has same number of columns as F matrix
%
% call lapack  DGEQRF() to perform QR factorization
%
% the Q matrix is encoded as Q = H(1)*H(2) ... H(k)
% so to get explicit version of Q,  
% we need to apply action of Q to columns of identity matrix
% by calling  lapack  DORMQR()
% ---------------------------------
[Q,R] = qr( F,0);
return;
end
