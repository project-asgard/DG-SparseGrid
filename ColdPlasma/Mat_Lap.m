% construct matrix for Lap 
global invM

Lap = GG-L-J+alpha*Q;

dofs = size(Lap,1);
II = speye(dofs);

Mat_Lap = kron(II,kron(II,Lap))+kron(II,kron(Lap,II))+kron(Lap,kron(II,II));

invMat = inv(Mat_Lap);

invM = blkdiag(invMat,invMat,invMat);