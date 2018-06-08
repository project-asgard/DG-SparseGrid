function [batch_list1, batch_list2, batch_list3, batch_list4, batch_list5, Y ] = ...
              kronmult5_batch(A1,A2,A3,A4,A5, X, ...
                  batch_list1, batch_list2, batch_list3, batch_list4, batch_list5)
% [batch_list1, batch_list2, batch_list3, batch_list4, batch_list5, Y ] = ...
%               kronmult5_batch(A1,A2,A3,A4,A5, X, ...
%                   batch_list1, batch_list2, batch_list3, batch_list4, batch_list5)

global idebug;

nrow1 = size(A1,1);
ncol1 = size(A1,2);

nrow2 = size(A2,1);
ncol2 = size(A2,2);

nrow3 = size(A3,1);
ncol3 = size(A3,2);

nrow4 = size(A4,1);
ncol4 = size(A4,2);

nrow5 = size(A5,1);
ncol5 = size(A5,2);

nrowX = ncol1*ncol2*ncol3*ncol4*ncol5;

isok = (mod(numel(X), nrowX) == 0);
if (~isok),
 error(sprintf('kronmult5_batch: numel(X)=%g, nrowX=%g', ...
                           numel(X),    nrowX ));
 return;
end;


nvec = numel( X )/nrowX;

nrowY = nrow1*nrow2*nrow3*nrow4*nrow5;


% -------------------------------------
% Ytmp = kronmult4( A2,A3,A4, A5, X );
% -------------------------------------
[batch_list1, batch_list2, batch_list3, batch_list4,  Ytmp] = ...
     kronmult4_batch( A2,A3,A4, A5, X, ...
         batch_list1, batch_list2, batch_list3, batch_list4 );

if (idebug >= 1),
  disp(sprintf('kronmult5_batch: numel(Ytmp)=%g', numel(Ytmp)));
end;

isok = (mod(numel(Ytmp),nvec) == 0);
if (~isok),
  error(sprintf('kronmult5_batch: numel(Ytmp)=%g, nvec=%g', ...
                            numel(Ytmp),    nvec ));
  return;
end;

% X is (ncol5*ncol4*ncol3*ncol2*ncol1) by nvec
% Ytmp = kronmult4( A2,A3,A4,A5,  X) so X appears to be (ncol5*ncol4*ncol3*ncol2) by (ncol1*nvec)
% Ytmp is (nrow5*nrow4*nrow3*nrow2) by (ncol1*nvec)
% so Ytmp can be reshaped as  (nrow5*nrow4*nrow3*nrow2) * ncol1  by nvec

nrowYtmp = (nrow5*nrow4*nrow3*nrow2) * ncol1;
Ytmp = reshape(Ytmp, [nrowYtmp, nvec]);

Y = zeros(nrowY, nvec);

for i=1:nvec,
 % -----------------------------------------------------------------
 % Yi = reshape( Ytmp(:,i), [nrowYtmp/ncol1, ncol1]) * transpose(A1);
 % -----------------------------------------------------------------
 mm = nrowYtmp/ncol1;
 kk = ncol1;
 nn = nrow1;
 transA = 'N';
 transB = 'T';
 alpha = 1;
 beta = 0;
 Amat = reshape( Ytmp(1:(mm*kk),i), mm, kk );
 Bmat = A1;
 Cmat = reshape( Y(1:(mm*nn),i), mm, nn );

  nbatch =  batch_list5.nbatch;
  nbatch = nbatch+1;
  batch_list5.nbatch = nbatch;

  batch_list5.mlist(nbatch) = mm;
  batch_list5.nlist(nbatch) = nn;
  batch_list5.klist(nbatch) = kk;

  batch_list5.transA(nbatch) = transA;
  batch_list5.transB(nbatch) = transB;

  batch_list5.alpha(nbatch) = alpha;
  batch_list5.beta(nbatch) = beta;

  batch_list5.Alist{nbatch} = Amat;
  batch_list5.Blist{nbatch} = Bmat;
  batch_list5.Clist{nbatch} = Cmat;

% ----------------------------------------
% perform computation to check correctness
% ----------------------------------------
  Cmat = gemm( transA, transB, mm, nn,kk, alpha, Amat, Bmat, beta, Cmat);

 Y(1:nrowY,i) = reshape( Cmat(1:mm,1:nn), [nrowY,1]);
end;

