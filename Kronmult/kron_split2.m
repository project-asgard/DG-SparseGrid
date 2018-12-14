function [iperm] = kron_split2( nA, nB1, nB2 )
%
% [iperm] = kron_split2( nA, nB1, nB2 )
% kron(A, [B1,B2]) * P = [ kron(A,B1) , kron(A,B2) ]
%
nB = nB1 + nB2;
ip = reshape( 1:(nB*nA), [nB,nA]);
%
% ----------------------------------------------
% ip(:,:) has indexing of the form "ip( ib, ia)"
% ----------------------------------------------
%
% ---------------------------------
% extract  the part with B1, and B2
% ---------------------------------
ip1 = ip(1:nB1, 1:nA);
ip2 = ip( (nB1+1):nB, 1:nA);

% -------------------------------
% generate new permutation vector
% -------------------------------
iperm = zeros( nB*nA,1);
iperm(1:(nB1*nA)) = reshape( ip1, (nB1*nA),1);
iperm( (nB1*nA) + (1:(nB2*nA)) ) = reshape( ip2, (nB2*nA),1);

