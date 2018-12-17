function [iperm] = kron_split( nA,  nBarray )
% -------------------------------------------------------
% [iperm] = kron_split( nA,  nBarray )
% 
% kron( A, [B1, ... Bn]) * P equals to  [kron(A,B1), ..., kron(A,Bn)]
% 
% similarly
% P* kron( A, [B1;
%              B2; 
%              ...
%              Bn]) 
% equals to 
%
% [ kron(A,B1);
%   kron(A,B2);
%   ...
%   kron(A,Bn) ]
% -------------------------------------------------------
idebug = 0;

numB = numel(nBarray);
nB = sum( nBarray(1:numB) );
if (idebug >= 1),
        disp(sprintf('kron_split: numB = %g, nB = %g ', numB, nB ));
        disp(sprintf('nBarray'));
        nBarray
end;


% ------------------------------------
% use recursion to call  kron_split2()
% ------------------------------------
if (numB == 1),
        % -------------------------------
        % trivial case, just the identity
        % -------------------------------
        iperm = reshape( 1:(nB*nA), (nB*nA),1);
        return;
end;

if (numB == 2),
        % -------------------------------------
        % only 2 sub blocks, use kron_split2()
        % -------------------------------------
        nB1 = nBarray(1);
        nB2 = nBarray(2);
        iperm = kron_split2( nA, nB1, nB2 );
        return;
end;


% ----------------------------------------------
% there are 3 or more sub blocks, need recursion
% ----------------------------------------------

% -----------------------------------------------------
% split   as nBarray(1:imid),  nBarray( (imid+1):numB )
% -----------------------------------------------------
imid = max(1,min(numB-1, floor(numB/2)));

iperm1 = kron_split( nA, nBarray( 1:imid) );
iperm2 = kron_split( nA, nBarray( (imid+1):numB ) );

nB1 = sum(nBarray(1:imid));
nB2 = sum(nBarray((imid+1):numB ));
iperm12 = kron_split2( nA, nB1, nB2 );

iperm = zeros( nA * (nB1 + nB2),1);

iptmp = iperm12( 1:(nA*nB1) );
iperm( 1:(nA*nB1) ) = iptmp( iperm1(1:(nA*nB1)) );


i1 = (nA*nB1) + 1;
i2 = (nA*nB1) + (nA*nB2);
iptmp = iperm12( i1:i2 );

iperm( i1:i2 ) = iptmp( iperm2(1:(nA*nB2)) );



















