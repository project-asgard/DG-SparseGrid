function [batch_list1, batch_list2, batch_list3, batch_list4, batch_list5, batch_list6, Y ] = ...
              kronmult6_batch(A1,A2,A3,A4,A5, A6, X, ...
                  batch_list1, batch_list2, batch_list3, batch_list4, batch_list5, batch_list6)
% [batch_list1, batch_list2, batch_list3, batch_list4, batch_list5, batch_list6, Y ] = ...
%               kronmult6_batch(A1,A2,A3,A4,A5, A6, X, ...
%                   batch_list1, batch_list2, batch_list3, batch_list4, batch_list5, batch_list6)

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

nrow6 = size(A6,1);
ncol6 = size(A6,2);

nrowX = ncol1*ncol2*ncol3*ncol4*ncol5*ncol6;

isok = (mod(numel(X), nrowX) == 0);
if (~isok),
 error(sprintf('kronmult6_batch: numel(X)=%g, nrowX=%g', ...
                           numel(X),    nrowX ));
 return;
end;


nvec = numel( X )/nrowX;
if (idebug >= 1),
 disp(sprintf('kronmult6_batch: (%d,%d) (%d,%d) (%d,%d) (%d,%d) (%d,%d) (%d,%d) nvec=%d', ...
          nrow1,ncol1, ...
          nrow2,ncol2, ...
          nrow3,ncol3, ...
          nrow4,ncol4, ...
          nrow5,ncol5, ...
          nrow6,ncol6, ...
          nvec ));
end;

nrowY = nrow1*nrow2*nrow3*nrow4*nrow5*nrow6;
Y = zeros(nrowY, nvec);
  
[flops1,flops2,imethod] = flops_kron6( nrow1,ncol1, nrow2,ncol2, ...
                                       nrow3,ncol3, nrow4,ncol4, ...
                                       nrow5,ncol5, nrow6,ncol6);
use_method_1  = (imethod == 1);
if (use_method_1),
  
  % -------------------------------------
  % Ytmp = kronmult5( A2,A3,A4, A5, A6, X );
  % -------------------------------------
  [batch_list1, batch_list2, batch_list3, batch_list4,  batch_list5,Ytmp] = ...
       kronmult5_batch( A2,A3,A4, A5, A6, X, ...
           batch_list1, batch_list2, batch_list3, batch_list4, batch_list5 );
  
  
  isok = (mod(numel(Ytmp),nvec) == 0);
  if (~isok),
    error(sprintf('kronmult6_batch: numel(Ytmp)=%g, nvec=%g', ...
                              numel(Ytmp),    nvec ));
    return;
  end;
  
  % -----------------------------------------------------
  % X is (ncol6*ncol5*ncol4*ncol3*ncol2*ncol1) by nvec
  % Ytmp = kronmult5( A2,A3,A4,A5,A6,   X) so X appears to be 
  % (ncol6*ncol5*ncol4*ncol3*ncol2) by (ncol1*nvec)
  % Ytmp is (nrow6*nrow5*nrow4*nrow3*nrow2) by (ncol1*nvec)
  % so Ytmp can be reshaped as  (nrow6*nrow5*nrow4*nrow3*nrow2) * ncol1  by nvec
  % -----------------------------------------------------
  
  nrowYtmp = (nrow6*nrow5*nrow4*nrow3*nrow2) * ncol1;
  Ytmp = reshape(Ytmp, [nrowYtmp, nvec]);
  
  
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
  
    nbatch =  batch_list6.nbatch;
    nbatch = nbatch+1;
    batch_list6.nbatch = nbatch;
  
    batch_list6.mlist(nbatch) = mm;
    batch_list6.nlist(nbatch) = nn;
    batch_list6.klist(nbatch) = kk;
  
    batch_list6.transA(nbatch) = transA;
    batch_list6.transB(nbatch) = transB;
  
    batch_list6.alpha(nbatch) = alpha;
    batch_list6.beta(nbatch) = beta;
  
    batch_list6.Alist{nbatch} = Amat;
    batch_list6.Blist{nbatch} = Bmat;
    batch_list6.Clist{nbatch} = Cmat;
  
  % ----------------------------------------
  % perform computation to check correctness
  % ----------------------------------------
    Cmat = gemm( transA, transB, mm, nn,kk, alpha, Amat, Bmat, beta, Cmat);
  
   Y(1:nrowY,i) = reshape( Cmat(1:mm,1:nn), [nrowY,1]);
  end; % end for

else

% -----------------------------------------------------
% Y = kron( A2, A3, A4, A5, A6 ) * (   X * transpose(A1) )
% -----------------------------------------------------
   X = reshape( X, nrowX, nvec );
   Ytmp = zeros( ncol2*ncol3*ncol4*ncol5*ncol6,   nrow1*nvec );
   

   for i=1:nvec,
     i1 = 1 + (i-1)*nrow1;
     i2 = i1 + nrow1 - 1;

%      ---------------------------
%      Xi = reshape( X(:,i),  (ncol2*ncol3*ncol4*ncol5*ncol6), ncol1 );
%      Ytmpi =  Xi(1:(ncol2*ncol3*ncol4*ncol5*ncol6),1:ncol1) * ...
%                   transpose( A1(1:nrow1,1:ncol1));
% 
% 
%      Ytmp(1:(ncol2*ncol3*ncol4*ncol5*ncol6),i1:i2) = ...
%            Ytmpi(1:(ncol2*ncol3*ncol4*ncol5*ncol6), 1:nrow1);
%      ---------------------------
       transA = 'N';
       transB = 'T';
       mm = (ncol2*ncol3*ncol4*ncol5*ncol6); 
       nn = nrow1;
       kk = ncol1;
       alpha = 1;
       beta = 0;

       Amat = reshape( X(:,i), mm,kk);
       Bmat = A1;
       Cmat = Ytmp(1:mm, i1:i2 );

       Cmat = gemm(transA, transB, mm,nn,kk, alpha, Amat, Bmat, beta, Cmat);
       Ytmp(1:mm, i1:i2) = Cmat(1:mm,1:nn);


       nbatch = batch_list1.nbatch + 1;
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

    end;


%     ----------------------------------
%     Y = kronmult5(A2,A3,A4,A5, Ytmp );
%     ----------------------------------
      [batch_list2, batch_list3, batch_list4, batch_list5, batch_list6, Y] = ...
                 kronmult5_batch( A2,A3,A4,A5,A6,  Ytmp, ...
                      batch_list2, batch_list3, batch_list4, batch_list5, batch_list6 );
             

end;


Y = reshape(Y, nrowY, nvec);
end
