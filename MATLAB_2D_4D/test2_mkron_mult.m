% 
rc = [2, 10; ...
      3, 11; ...
      4, 12; ...
      5, 13; ...
      15, 7 ];
rc = transpose(rc);

nkron = size(rc,2);
for i=1:nkron,
  nrow = rc(1,i);
  ncol = rc(2,i);
  Acell{i} = rand( nrow, ncol );
end;

nvec = 3;
X = rand( prod(rc(2,1:nkron)), nvec );

t1 = tic;
Y1 = kron_multd( nkron, Acell, X );
disp(sprintf('1st kron_mult took %g sec', toc(t1)));

t1 = tic;
Y1 = kron_multd( nkron, Acell, X );
disp(sprintf('2nd kron_mult_took %g sec', toc(t1)));

t1 = tic;
Y2 = mkron_multd( nkron, Acell, X );
disp(sprintf('1st mkron_mult took %g sec', toc(t1)));

t1 = tic;
Y2 = mkron_multd( nkron, Acell, X );
disp(sprintf('2nd mkron_mult_took %g sec', toc(t1)));

err = abs(Y1-Y2);
err = max( abs(err(:)) );
disp(sprintf('max error is %g ', err ));
