% ----------------------------------
% simple test of kron_mult routines
% ----------------------------------
nvec = 3;
Acell{1} = rand(2,3);
Acell{2} = rand(3,4);
Acell{3} = rand(4,5);
Acell{4} = rand(6,7);
Acell{5} = rand(5,8);
Acell{6} = rand(2,2);
nkron = 6;

rc = zeros(2,nkron);
for k=1:nkron,
  rc(1,k) = size( Acell{k},1);
  rc(2,k) = size( Acell{k},2);
end;

nkron = 6;
for k=2:nkron,
  isizeX = prod(rc(2,1:k));
  X = rand( isizeX,nvec);
  if (k==2),
    Y = kron( Acell{1},Acell{2} ) * X;
    Ykron = kron_mult2( Acell{1},Acell{2},X);
  elseif (k==3),
    Y = kron( Acell{1},kron(Acell{2},Acell{3}) ) * X;
    Ykron = kron_mult3( Acell{1},Acell{2},Acell{3},X);
  elseif (k==4),
    A12 = kron( Acell{1},Acell{2} );
    A34 = kron( Acell{3},Acell{4} );
    Y = kron( A12, A34 ) * X;
    Ykron = kron_mult4( Acell{1},Acell{2},Acell{3},Acell{4},X);
  elseif (k==5),
    A12 = kron( Acell{1},Acell{2} );
    A34 = kron( Acell{3},Acell{4} );
    Y = kron( kron(A12, A34), Acell{5} ) * X;
    Ykron = kron_mult5( Acell{1},Acell{2},Acell{3},Acell{4},Acell{5},X);
  elseif (k == 6),
    A12 = kron( Acell{1},Acell{2} );
    A34 = kron( Acell{3},Acell{4} );
    A56 = kron( Acell{5},Acell{6} );
    Y = kron( kron(A12,A34), A56 ) * X;
    Ykron = kron_mult6( Acell{1},Acell{2},Acell{3},Acell{4},Acell{5},Acell{6},X);
  end;

    err = norm(Y(:) - Ykron(:));
    if (err > 1.0e-6),
      disp(sprintf('k=%g, err=%g',k,err));
    end;
end;




