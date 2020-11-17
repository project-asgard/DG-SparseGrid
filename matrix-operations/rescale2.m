function [As,D0] = rescale2( A )
% As = rescale2(A)
% row and column scaling so that
% after rescaling, diagonal entries
% are about unit value
% -----------------------------
dd = sqrt( abs(diag(A)) );
dd = dd + sign(dd)*eps;
D0=spdiag( 1./dd, 0);
As = D0 * A * D0;
