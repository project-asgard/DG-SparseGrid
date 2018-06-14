global idebug;
idebug = 1;

A1 = rand(2,3);
A2 = rand(3,4);
A3 = rand(4,5);
A4 = rand(5,6);
A5 = rand(6,7);

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

nvec = 2;
tol = 1e-7;
nerr = 0;
% -------------------------
% check   Y = kron(A1,A2)*X
% -------------------------
X2  = rand( ncol1*ncol2, nvec);
Y2ok = kron(A1,A2)*X2;
Y2 = kronmult2(A1,A2,X2);
diff = max(abs(Y2ok(:)-Y2(:)));
isok = (diff < tol * numel(Y2ok) );
if (~isok),
  disp(sprintf('kronmult2 failed, diff=%g',diff));
  nerr = nerr + 1;
end;

% --------------------------
% check Y = kron(A1,A2,A3)*X
% --------------------------
X3  = rand( ncol1*ncol2*ncol3, nvec);
Y3ok = kron(A1,kron(A2,A3))*X3;
Y3 = kronmult3(A1,A2,A3,X3);
diff = max(abs(Y3ok(:)-Y3(:)));
isok = (diff < tol * numel(Y3ok) );
if (~isok),
  disp(sprintf('kronmult3 failed, diff=%g',diff));
  nerr = nerr + 1;
end;


% --------------------------
% check Y = kron(A1,A2,A3,A4)*X
% --------------------------
X4  = rand( ncol1*ncol2*ncol3*ncol4, nvec);
Y4ok = kron(A1,kron(A2,kron(A3,A4)))*X4;
Y4 = kronmult4(A1,A2,A3,A4,  X4);
diff = max(abs(Y4ok(:)-Y4(:)));
isok = (diff < tol * numel(Y4ok) );
if (~isok),
  disp(sprintf('kronmult4 failed, diff=%g',diff));
  nerr = nerr + 1;
end;

% --------------------------
% check Y = kron(A1,A2,A3,A4,A5)*X
% --------------------------
X5  = rand( ncol1*ncol2*ncol3*ncol4*ncol5, nvec);
Y5ok = kron(A1,kron(A2,kron(A3,kron(A4,A5))))*X5;
Y5 = kronmult5(A1,A2,A3,A4,A5,  X5);
diff = max(abs(Y5ok(:)-Y5(:)));
isok = (diff < tol * numel(Y5ok) );
if (~isok),
  disp(sprintf('kronmult5 failed, diff=%g',diff));
  nerr = nerr + 1;
end;

if (nerr == 0),
 disp('All OK');
end;





