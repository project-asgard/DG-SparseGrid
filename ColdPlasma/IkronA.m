function B = IkronA(A,nz)
% nz is the size of eye matrix


n1 = size(A,1);
B = sparse(n1*nz);

[i,j,aij] = find(A);
nAx = size(i,1);

indexJ = ones(nAx,1)*([0:nz-1]*nz);

II = i*ones(1,nz)+indexJ;
JJ = j*ones(1,nz)+indexJ;
Aij = aij*ones(1,nz);

B = sparse(II,JJ,Aij);

