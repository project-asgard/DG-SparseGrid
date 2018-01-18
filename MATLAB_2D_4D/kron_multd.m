function Y = kron_multd( nkron, Acell, X )
%  Y = kron_multd( Acell, X )
%  generic multi-dimension kron product
%
% reference implementation
%
idebug = 0;
isok = (length(Acell) >= nkron);
if (~isok),
  error(sprintf('kron_multd: invalid nkron=%d, length(Acell)=%d', ...
        nkron, length(Acell) ));
  return;
end;

  rc = zeros(2,nkron);
  for k=1:nkron,
    rc(1,k) = size( Acell{k},1);
    rc(2,k) = size( Acell{k},2);
  end;
  
  isizeX = prod( rc(2,1:nkron) );
  isizeY = prod( rc(1,1:nkron) );
  
  nvec = prod(size(X))/isizeX;
  isok = (mod( prod(size(X)), isizeX) == 0);
  if (~isok),
    error(sprintf('kron_multd:nkron=%d,prod(size(X))=%g,isizeX=%g', ...
               nkron, prod(size(X)), isizeX ));
  end;

if (nkron == 1),
  nc = size(Acell{1},2);
  Y = Acell{1} * reshape(X, [nc, prod(size(X))/nc]);
else
  % ------------------------------
  % general case require recursion
  % ------------------------------
  use_split_first = 1;

  if (use_split_first),
  % -------------------------------
  % kron(A(1), ..., A(nkron)) * X
  %
  % use simplified fixed algorithm
  % computed as
  %
  % kron( A(2)... A(n-1))*X * transpose(A(1))
  % -------------------------------
  n = nkron;
  nm1 = n-1;

  m1 = rc(1,1);
  n1 = rc(2,1);
  m2 = prod( rc(1,2:n));
  n2 = prod( rc(2,2:n));

  for i=1:nm1,
    Acell_tmp{i} = Acell{1+i};
  end;

  X = reshape( X, n2, n1*nvec );
  Z = kron_multd( nm1, Acell_tmp, X ); % note Z is m2 by (n1*nvec)
  Z = reshape( Z, m2*n1, nvec );
  Y = zeros( m2*m1, nvec );
  A1 = Acell{1}; % note An is m1 by n1

  if (idebug >= 1),
    disp(sprintf('kron_multd: nkron=%d,m1=%d,n1=%d,m2=%d,n2=%d', ...
                              nkron,   m1,   n1,   m2,   n2 ));
    disp(sprintf('size(A1)=(%d,%d), size(Z)=(%d,%d)', ...
                  size(A1,1),size(A1,2),  size(Z,1),size(Z,2) ));
  end;


  for i=1:nvec,
    Zi = reshape( Z(:,i), m2,n1);
    Yi = Zi * transpose( A1 ); % note Yi is m2 by m1
    Y(:,i) = reshape( Yi, m2*m1,1);
  end;

  else  
  
  % -------------------------------
  % kron(A(1), ..., A(nkron)) * X
  %
  % use simplified fixed algorithm
  % computed as
  %
  % kron( kron(A(1)..A(nkron-1)),  A(nkron))*X
  % A(nkron)*X*transpose( kron(A(1),...A(nkron-1))
  % or
  % Z = A(nkron)*X, then
  %
  % Z * transpose( kron(A(1)...,A(nkron-1))
  %
  % transpose( (  kron(A(1)...A(nkron-1))*transpose(Z)  )
  % -------------------------------
  X = reshape( X, isizeX,nvec);


  for k=1:nvec,
  
   Xk = X(:,k);

   nc = size(Acell{nkron},2);
   Z = Acell{nkron} * reshape( Xk, [nc, prod(size(Xk))/nc]);
 
   nc2 = prod( rc(2,1:(nkron-1)) );
   isok = mod(prod(size(Z)),nc2) == 0;
   if (~isok),
     error(sprintf('kron_multd: invalid size(Z)=%g, size(Acell{nkron},2)=%g', ...
           prod(size(Z)), size(Acell{nkron},2) ));
    return;
   end;

   Ytmp = kron_multd( nkron-1,Acell, transpose(Z));

  
   Y(:,k) = reshape(transpose(Ytmp), [isizeY,1]);
  end;

 end;
end;
