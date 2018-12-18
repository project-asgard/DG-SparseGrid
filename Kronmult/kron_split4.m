
function [iperm] = kron_split4( nA, nB1, nB2, nB3, nB4 )
% ---------------------------------------------------
% [iperm] = kron_split4( nA, nB1, nB2, nB3, nB4 )
% kron(A, [B1,B2,B3,B4]) * P equals [kron(A,B1),kron(A,B2),kron(A,B3),kron(A,B4)]
%
% similarly
% P * kron(A, [B1;
%              B2;
%              B3;
%              B4])
% equals
% [ kron(A,B1); 
%   kron(A,B2);
%   kron(A,B3);
%   kron(A,B4) ]
% ---------------------------------------------------

iperm = kron_split( nA, [nB1, nB2, nB3, nB4] );
