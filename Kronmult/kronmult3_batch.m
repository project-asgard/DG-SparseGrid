function [batch_list1,batch_list2,batch_list3,Y] = ...
        kronmult3_batch(A1,A2,A3, X, ...
                        batch_list1,batch_list2,batch_list3)
% [batch_list1,batch_list2,batch_list3,Y] = ...
%        kronmult3_batch(A1,A2,A3, X, ...
%                        batch_list1,batch_list2,batch_list3)
global idebug;

nrow1 = size(A1,1);
ncol1 = size(A1,2);

nrow2 = size(A2,1);
ncol2 = size(A2,2);

nrow3 = size(A3,1);
ncol3 = size(A3,2);

nrowX = ncol1*ncol2*ncol3;

isok = (mod( numel(X), nrowX ) == 0);
if (~isok),
 error(sprintf('kronmult3_batch: numel(X)=%g, nrowX=%g', ...
                                numel(X),    nrowX));
 return;
end;

nvec = numel(X)/nrowX;

if (idebug >= 1),
 disp(sprintf('kronmult3_batch: (%d,%d) (%d,%d) (%d,%d), nvec=%d',...
        nrow1,   ncol1,    nrow2,   ncol2,    nrow3,   ncol3,    nvec));
end;

nrowY = nrow1*nrow2*nrow3;
Y = zeros(nrowY, nvec);

[flops1,flops2,imethod] = flops_kron3( nrow1,ncol1, nrow2,ncol2, nrow3,ncol3);
use_method_1  = (imethod == 1);
if (use_method_1),

  % --------------------------
  % Ytmp = kronmult2( A2,A3,X);
  % --------------------------
  [batch_list1, batch_list2, Ytmp] = kronmult2_batch( A2, A3, X, ...
                                           batch_list1, batch_list2);
  
  isok = (mod( numel(Ytmp), nvec ) == 0);
  if (~isok),
    error(sprintf('kronmult3_batch: numel(Ytmp)=%g, nvec=%g', ...
                              numel(Ytmp),    nvec));
    return;
  end;
  
  % -----------------------------------------------
  % X is ncol1*ncol2*ncol3 by nvec
  %
  % Ytmp = kronmult2(A2,A3,  X), X appears as (ncol3*ncol2) by (ncol1*nvec)
  % so Ytmp is  (nrow2*nrow3) by (ncol1*nvec)
  % and Ytmp can be reshaped to be (nrow2*nrow3*ncol1) by nvec
  % -----------------------------------------------
  %nrowYtmp = numel(Ytmp)/nvec;
  nrowYtmp = (nrow2*nrow3*ncol1);
  
  Ytmp = reshape( Ytmp, [nrowYtmp, nvec] );
  
  
  % -----------------------------------------------------
  % note: may be task parallelism or batch gemm operation
  % -----------------------------------------------------
  isok = (mod(nrowYtmp,ncol1) == 0);
  if (~isok),
    error(sprintf('kronmult3: nrowYtmp=%g, ncol1=%g', ...
                              nrowYtmp,    ncol1 ));
    return;
  end;
  msize  = nrowYtmp/ncol1;
  
  for i=1:nvec,
    % ------------------------------------------------------
    % Yi = reshape(Ytmp(:,i), [msize, ncol1])*transpose(A1);
    % ------------------------------------------------------
  
    mm = msize;
    kk = ncol1;
    nn = nrow1;
    transA = 'N';
    transB = 'T';
    alpha = 1;
    beta = 0;
    
    Amat = reshape( Ytmp(1:(mm*kk),i), mm, kk);
    Bmat = A1;
    Cmat = reshape( Y(1:(mm*nn),i), mm,nn );
  
  
    nbatch = batch_list3.nbatch;
    nbatch = nbatch+1;
    batch_list3.nbatch = nbatch;
    
    batch_list3.mlist(nbatch) = mm;
    batch_list3.nlist(nbatch) = nn;
    batch_list3.klist(nbatch) = kk;
  
    batch_list3.transA(nbatch) = transA;
    batch_list3.transB(nbatch) = transB;
  
    batch_list3.alpha(nbatch) = alpha;
    batch_list3.beta(nbatch) = beta;
  
    batch_list3.Alist{nbatch} = Amat;
    batch_list3.Blist{nbatch} = Bmat;
    batch_list3.Clist{nbatch} = Cmat;
  
  
  % ----------------------------------------
  % perform computation to check correctness
  % ----------------------------------------
    Cmat = gemm( transA, transB, mm, nn,kk, alpha, Amat, Bmat, beta, Cmat);
  
    Y(1:nrowY,i) = reshape(Cmat(1:mm,1:nn), nrowY,1);
  end; % end for


else

%  ---------------------------------------
%  Y = kron( A2, A3) ( X * transpose(A1) )
%  ---------------------------------------
%    Xi = zeros( ncol2*ncol3, ncol1 );
%    Ytmpi = zeros( ncol2*ncol3, nrow1 );
   X = reshape( X, nrowX, nvec );
   Ytmp = zeros( (ncol2*ncol3), nrow1*nvec );
   for i=1:nvec,
     % --------------------------------
     % note  nrowX = ncol1*ncol2*ncol3;
     % so nrowX/ncol1 =  ncol2 * ncol3
     % --------------------------------
     i1 = 1 + (i-1)*nrow1;
     i2 = i1 + nrow1 - 1;

%      -----------------------------------------
%      Xi = reshape( X(:,i), (ncol2 * ncol3), ncol1 );
%      Ytmpi(1:(ncol2*ncol3), 1:nrow1) = ...
%        Xi(1:(ncol2*ncol3), 1:ncol1) * transpose(A1(1:nrow1,1:ncol1));
%      Ytmp(:,i1:i2) = Ytmpi(1:(ncol2*ncol3), 1:nrow1 );
%      -----------------------------------------
       transA = 'N';
       transB = 'T';
       mm = ncol2*ncol3;
       nn = nrow1;
       kk = ncol1;
       alpha = 1;
       beta = 0;
      
       Amat = reshape( X(:,i), (ncol2*ncol3), ncol1 );
       Bmat = A1;
       Cmat  = Ytmp(1:(ncol2*ncol3), i1:i2);

       Cmat = gemm( transA, transB, mm,nn,kk, alpha, Amat, Bmat, beta, Cmat);
       Ytmp(1:(ncol2*ncol3), i1:i2) = Cmat(1:mm,1:nn);


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

   end; % end for 

%    ------------------------------
%    Y = kronmult2( A2, A3, Ytmp );
%    ------------------------------
   [batch_list2, batch_list3, Y] = kronmult2_batch( A2, A3, Ytmp, batch_list2, batch_list3 );
end;

Y = reshape(Y, nrowY, nvec );
end
