% test kron_split
A = rand(3,4);
nrowA = size(A,1);
ncolA = size(A,2);

B1 = rand(5,5);
B2 = rand(5,6);
B3 = rand(5,7);
B4 = rand(5,8);
B5 = rand(5,9);

B = [B1,B2,B3,B4,B5];
AxB = kron(A,B);

AxB1 = kron(A,B1);
AxB2 = kron(A,B2);
AxB3 = kron(A,B3);
AxB4 = kron(A,B4);
AxB5 = kron(A,B5);

nBarray = [ size(B1,2), size(B2,2), size(B3,2), size(B4,2), size(B5,2)];

iperm_col = kron_split( ncolA,  nBarray );

diff_col = norm( AxB(:,iperm_col) - [AxB1,AxB2,AxB3,AxB4,AxB5],1);
disp(sprintf('diff of kron_split on columns = %g ', diff_col ));

C1 = rand(5,5);
C2 = rand(6,5);
C3 = rand(7,5);
C4 = rand(8,5);
C5 = rand(9,5);

C = [C1; ...
     C2; ...
     C3; ...
     C4; ...
     C5];
AxC = kron(A,C);

AxC1 = kron(A,C1);
AxC2 = kron(A,C2);
AxC3 = kron(A,C3);
AxC4 = kron(A,C4);
AxC5 = kron(A,C5);

nCarray = [size(C1,1), size(C2,1), size(C3,1), size(C4,1), size(C5,1) ];
iperm_row = kron_split( nrowA, nCarray );

diff_row = norm( AxC( iperm_row,:) - [ AxC1; ...
                                       AxC2; ...
                                       AxC3; ...
                                       AxC4; ...
                                       AxC5 ], 1 );
disp(sprintf('diff of kron_split on rows = %g ', diff_row ));
