function [batch_list1,Y] = kronmult1_batch(A1, X, batch_list1 )
% [batch_list1,Y] = kronmult1_batch(A1, X, batch_list1 )
% generate  list for dgemm_vbatch
%
global idebug;

nrow1 = size(A1,1);
ncol1 = size(A1,2);


% -----------
% extra check
% -----------
isok = (mod( numel(X), ncol1) == 0);
if (~isok),
  error(sprintf('kronmult1_batch: numel(X)=%g, ncol1=%g', ...
                         numel(X), ncol1 ));
  return;
end;

nvec = numel(X)/ncol1;

if (idebug >= 1),
  disp(sprintf('kronmult1_batch: (%d,%d) nvec=%d', ...
                    nrow1,    ncol1,    nvec));
end;

Y = zeros( nrow1, nvec );

% --------------------------------------------
% Y = A1 * reshape(X, ncol1, numel(X)/ncol1 );
% --------------------------------------------
transA = 'N';
transB = 'N';
mm = nrow1;
kk = ncol1;
nn = nvec;

Amat = A1;
Bmat = reshape(X,kk,nn);
alpha = 1;
beta = 0;
Cmat = Y(1:mm,1:nn);

nbatch = batch_list1.nbatch;
nbatch = nbatch+1;
batch_list1.nbatch = nbatch;

batch_list1.mlist(nbatch) = mm;
batch_list1.nlist(nbatch) = nn;
batch_list1.klist(nbatch) = kk;

batch_list1.transA(nbatch) = transA;
batch_list1.transB(nbatch) = transB;

batch_list1.alpha(nbatch) = alpha;
batch_list1.beta(nbatch) = beta;

batch_list1.Alist{nbatch} = Amat;
batch_list1.Blist{nbatch} = Bmat;
batch_list1.Clist{nbatch} = Cmat;

% ----------------------------------------
% perform computation to check correctness
% ----------------------------------------
Cmat = gemm( transA, transB, mm, nn,kk, alpha, Amat, Bmat, beta, Cmat);
Y(1:mm,1:nn) = Cmat(1:mm,1:nn);



end
