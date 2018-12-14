A = rand(3,4);
B1 = rand(5,6);
B2 = rand(5,7);

nA  = size(A,2);
nB1 = size(B1,2);
nB2 = size(B2,2);

AxB  = kron( A, [B1,B2]);
AxB1 = kron(A, B1);
AxB2 = kron(A, B2);

iperm_col = kron_split2( nA, nB1, nB2 );

% ------------------------------------------------------------
% check   kron(A, [B1, B2]) * P ==  [ kron(A,B1), kron(A,B2) ]
% ------------------------------------------------------------
diff_col = norm(AxB(:,iperm_col) - [AxB1, AxB2],1);
disp(sprintf('diff of kron(A,[B1,B2])*P - [kron(A,B1),kron(A,B2)] = %g ', diff_col ));

% ----------------------------
% similarly test row partition
% ----------------------------
C1 = rand(6,5);
C2 = rand(7,5);
C = [C1; ...
     C2];

AxC  = kron(A,C);
AxC1 = kron(A,C1);
AxC2 = kron(A,C2);

nrowA  = size(A,1);
nrowC1 = size(C1,1);
nrowC2 = size(C2,1);

iperm_row = kron_split2( nrowA, nrowC1, nrowC2 );
diff_row = norm( AxC(iperm_row,:) - [ AxC1; ...
                                      AxC2], 1);
disp(sprintf('diff of P*kron(A,[C1;C2]) - [kron(A,C1); kron(A,C2)] = %g', diff_row));


