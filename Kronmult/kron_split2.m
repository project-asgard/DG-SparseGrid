function [iperm] = kron_split2( nA, nB1, nB2 )
% ---------------------------------------------------
% [iperm] = kron_split2( nA, nB1, nB2 )
% kron(A, [B1,B2]) * P equals [ kron(A,B1) , kron(A,B2) ]
%
% similarly
% P * kron(A, [B1;
%              B2])
% equals
% [ kron(A,B1); 
%   kron(A,B2)]
% ---------------------------------------------------
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
iperm(1:(nB1*nA)) = reshape( ip1, numel(ip1),1);

i1 = (nB1*nA) + 1;
i2 = (nB1*nA) + (nB2*nA);
iperm( i1:i2 ) = reshape( ip2, numel(ip2),1);

