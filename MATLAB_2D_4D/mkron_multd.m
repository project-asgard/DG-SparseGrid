function Y = mkron_multd( nkron, Acell, X )
%  Y = mkron_multd( nkron, Acell, X )
%
%  generic multi-dimension kron product
%
%  implementation to minimize flops
%
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
  [flops,isplit,imethod] = kron_minflops( rc ); 
  k = isplit;
  kp1 = k+1;
  n = nkron;
  n1 = prod(rc(2,1:k));
  n2 = prod(rc(2,kp1:n));
  m1 = prod(rc(1,1:k));
  m2 = prod(rc(1,kp1:n));

  X = reshape( X, n1*n2, nvec);
  Y = zeros( m1*m2,nvec );
  if (imethod == 1),
    %  X = reshape(X, n2,n1)
    %  T1 = kron( Akp1...An) *  X
    %  T1 is m2 by n1
    %  T1t = reshape( transpose(T1), n1, m2 )
    %  Y = transpose(kron( A1..Ak) * T1t)
  
    
    for ivec=1:nvec,
       Xi = reshape( X(:,ivec), n2,n1);

       % --------------------------
       % copy Acell1 = Acell{kp1:n}
       % --------------------------
       for j=1:(n-kp1+1),
          Acell1{j} = Acell{k+j};
       end;
       T1 = mkron_multd( (n-kp1+1), Acell1, Xi);

       T1t = reshape( transpose(T1), n1,m2);
       Yit = mkron_multd( k, Acell, T1t); % Y1t is m1 by m2
       Yi =  reshape(transpose(Yit), m2,m1); 
       Y(1:(m1*m2),ivec) = reshape(Yi, m1*m2,1);
    end;
  else
    % Xt = reshape( transpose(X), n1,n2);
    % T2 = kron( A1..Ak ) * Xt, T2 is m1 by n2
    % T2t = reshape(transpose(T2), n2, m1 )
    % Y = kron( Akp1..An) * T2t
   
    for ivec=1:nvec,
      Xi = reshape(X(:,ivec), n2,n1);
      Xit = reshape(transpose(Xi), n1, n2);

      T2 = mkron_multd( k, Acell, Xit );
      T2t = reshape( transpose(T2), n2, m1 );
      
      % --------------------------
      % copy Acell2 = Acell{kp1:n}
      % --------------------------
      for j=1:(n-kp1+1),
         Acell2{j} = Acell{k+j};
      end;
      Yi = mkron_multd( (n-kp1+1), Acell2, T2t ); % Yi is m2 by m1

      Y(1:(m1*m2),ivec) = reshape( Yi, m1*m2,1);
    end;
  end; % if (imethod)
 end; % if (nkron == 1)
end;
