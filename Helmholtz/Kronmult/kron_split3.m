
function [iperm] = kron_split3( nA, nB1, nB2, nB3 )
% ---------------------------------------------------
% [iperm] = kron_split3( nA, nB1, nB2, nB3 )
% kron(A, [B1,B2,B3]) * P equals [ kron(A,B1) , kron(A,B2), kron(A,B3) ]
%
% similarly
% P * kron(A, [B1;
%              B2;
%              B3])
% equals
% [ kron(A,B1); 
%   kron(A,B2);
%   kron(A,B3)]
% ---------------------------------------------------

iperm = kron_split( nA, [nB1, nB2, nB3] );
