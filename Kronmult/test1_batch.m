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

batch_list1 = clear_batch_list();
batch_list2 = clear_batch_list();
batch_list3 = clear_batch_list();
batch_list4 = clear_batch_list();
batch_list5 = clear_batch_list();

% --------------
% check Y1 = A1*X1
% --------------
disp('========================');
X1  = rand(ncol1,nvec);
Y1ok = A1*X1;
[batch_list1, Y1] = ...
       kronmult1_batch(A1, X1, ...
           batch_list1 );

diff = max(abs(Y1ok(:)-Y1(:)));
isok = (diff < tol * numel(Y1ok) );
if (~isok),
  disp(sprintf('kronmult1_batch failed, diff=%g',diff));
  nerr = nerr + 1;
end;


batch_list1 = clear_batch_list();
batch_list2 = clear_batch_list();
batch_list3 = clear_batch_list();
batch_list4 = clear_batch_list();
batch_list5 = clear_batch_list();
% -------------------------
% check   Y = kron(A1,A2)*X
% -------------------------
disp('========================');
X2  = rand( ncol1*ncol2, nvec);
Y2ok = kron(A1,A2)*X2;
[batch_list1, batch_list2, Y2 ] = ...
    kronmult2_batch(A1,A2,X2, ...
       batch_list1, batch_list2);
diff = max(abs(Y2ok(:)-Y2(:)));
isok = (diff < tol * numel(Y2ok) );
if (~isok),
  disp(sprintf('kronmult2_batch failed, diff=%g',diff));
  nerr = nerr + 1;
end;



batch_list1 = clear_batch_list();
batch_list2 = clear_batch_list();
batch_list3 = clear_batch_list();
batch_list4 = clear_batch_list();
batch_list5 = clear_batch_list();
% --------------------------
% check Y = kron(A1,A2,A3)*X
% --------------------------
disp('========================');
X3  = rand( ncol1*ncol2*ncol3, nvec);
Y3ok = kron(A1,kron(A2,A3))*X3;
[batch_list1, batch_list2, batch_list3,Y3] = ...
   kronmult3_batch(A1,A2,A3,X3, ...
      batch_list1, batch_list2, batch_list3);

diff = max(abs(Y3ok(:)-Y3(:)));
isok = (diff < tol * numel(Y3ok) );
if (~isok),
  disp(sprintf('kronmult3_batch failed, diff=%g',diff));
  nerr = nerr + 1;
end;


batch_list1 = clear_batch_list();
batch_list2 = clear_batch_list();
batch_list3 = clear_batch_list();
batch_list4 = clear_batch_list();
batch_list5 = clear_batch_list();
% --------------------------
% check Y = kron(A1,A2,A3,A4)*X
% --------------------------
disp('========================');
X4  = rand( ncol1*ncol2*ncol3*ncol4, nvec);
Y4ok = kron(A1,kron(A2,kron(A3,A4)))*X4;
[batch_list1, batch_list2, batch_list3, batch_list4, Y4] = ...
   kronmult4_batch(A1,A2,A3,A4,  X4, ...
      batch_list1, batch_list2, batch_list3, batch_list4);

diff = max(abs(Y4ok(:)-Y4(:)));
isok = (diff < tol * numel(Y4ok) );
if (~isok),
  disp(sprintf('kronmult4_batch failed, diff=%g',diff));
  nerr = nerr + 1;
end;


batch_list1 = clear_batch_list();
batch_list2 = clear_batch_list();
batch_list3 = clear_batch_list();
batch_list4 = clear_batch_list();
batch_list5 = clear_batch_list();
% --------------------------
% check Y = kron(A1,A2,A3,A4,A5)*X
% --------------------------
disp('========================');
X5  = rand( ncol1*ncol2*ncol3*ncol4*ncol5, nvec);
Y5ok = kron(A1,kron(A2,kron(A3,kron(A4,A5))))*X5;
[batch_list1, batch_list2, batch_list3, batch_list4, batch_list5, Y5]  = ...
    kronmult5_batch(A1,A2,A3,A4,A5,  X5, ...
    batch_list1, batch_list2, batch_list3, batch_list4, batch_list5);
diff = max(abs(Y5ok(:)-Y5(:)));
isok = (diff < tol * numel(Y5ok) );
if (~isok),
  disp(sprintf('kronmult5_batch failed, diff=%g',diff));
  nerr = nerr + 1;
end;

if (nerr == 0),
 disp('All OK');
end;

[total_flops1, max_flops1, min_flops1]  = flops_batch_list( batch_list1 );
[total_flops2, max_flops2, min_flops2]  = flops_batch_list( batch_list2 );
[total_flops3, max_flops3, min_flops3]  = flops_batch_list( batch_list3 );
[total_flops4, max_flops4, min_flops4]  = flops_batch_list( batch_list4 );
[total_flops5, max_flops5, min_flops5]  = flops_batch_list( batch_list5 );

nbatch1 = batch_list1.nbatch;
nbatch2 = batch_list2.nbatch;
nbatch3 = batch_list3.nbatch;
nbatch4 = batch_list4.nbatch;
nbatch5 = batch_list5.nbatch;

disp(sprintf('nbatch1=%d, total_flops1=%e, max_flops1=%e, min_flops=%e', ...
              nbatch1,    total_flops1,    max_flops1,    min_flops1 ));

disp(sprintf('nbatch2=%d, total_flops2=%e, max_flops2=%e, min_flops=%e', ...
              nbatch2,    total_flops2,    max_flops2,    min_flops2 ));

disp(sprintf('nbatch3=%d, total_flops3=%e, max_flops3=%e, min_flops=%e', ...
              nbatch3,    total_flops3,    max_flops3,    min_flops3 ));

disp(sprintf('nbatch4=%d, total_flops4=%e, max_flops4=%e, min_flops=%e', ...
              nbatch4,    total_flops4,    max_flops4,    min_flops4 ));


disp(sprintf('nbatch5=%d, total_flops5=%e, max_flops5=%e, min_flops=%e', ...
              nbatch5,    total_flops5,    max_flops5,    min_flops5 ));


disp(sprintf('total work = %e ', ...
    total_flops1 + total_flops2 + total_flops3 + total_flops4 + total_flops5 ));

