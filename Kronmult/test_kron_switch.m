
nrA = 3;
ncA = 4;
nrB = 5;
ncB = 6;
A = rand(nrA,ncA);
B = rand(nrB,ncB);

AxB = kron(A,B);
BxA = kron(B,A);
[ip_r,ip_c] = kron_switch( nrA,ncA, nrB,ncB);

diff =  norm( AxB( ip_r, ip_c ) - BxA,1);
disp(sprintf('norm(AxB(ip_r,ip_c)-BxA) = %g', diff));
