function B = AkronI(A,nz)
% nz is the size of eye matrix

[i,j,aij] = find(A);
II = (i-1)*nz+[1:nz];
JJ = (j-1)*nz+[1:nz];
Aij = aij*ones(1,nz);

B = sparse(II,JJ,Aij);