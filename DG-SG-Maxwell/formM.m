function M = formM( B )
% M = formM( B )
% --------------------
% M = [ 0,    M12,  M13;
%      -M12,  0,    M23;
%      -M13, -M23,  0]
%
% ---------------------
% M12 = kron(eye,eye,B)
% M13 = -kron(eye,B,eye)
% M23 = kron(B,eye,eye)
% ---------------------
n = size(B,1);
E = speye( n,n);
Z = sparse(n^3, n^3 );
M12 =  kron(E,kron(E,B));
M13 = -kron(E,kron(B,E));
M23 =  kron(B,kron(E,E));

M = [ Z,    M12, M13; ...
     -M12,  Z,   M23; ...
     -M13, -M23, Z ];
