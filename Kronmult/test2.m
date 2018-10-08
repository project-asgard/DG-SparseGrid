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


size_X = (ncol1 * ncol2 * ncol3 * ncol4 * ncol5) * nvec;
size_Y = (nrow1 * nrow2 * nrow3 * nrow4 * nrow5) * nvec;

X = rand( ncol1*ncol2*ncol3*ncol4*ncol5, nvec );
Y = zeros( nrow1 * nrow2 * nrow3 * nrow4 * nrow5, nvec );

disp(sprintf('nrow1=%g, ncol1=%g', nrow1,ncol1 ));
disp(sprintf('nrow2=%g, ncol2=%g', nrow2,ncol2 ));
disp(sprintf('nrow3=%g, ncol3=%g', nrow3,ncol3 ));
disp(sprintf('nrow4=%g, ncol4=%g', nrow4,ncol4 ));
disp(sprintf('nrow5=%g, ncol5=%g', nrow5,ncol5 ));
disp(sprintf('nvec=%g', nvec ));

disp(sprintf('size_X=%g, size_Y=%g ',  ...
              size_X,    size_Y ));


[total_mem_use, mem_use5] = mem_kron5( ...
                              nrow1, ncol1, ...
                              nrow2, ncol2, ...
                              nrow3, ncol3, ...
                              nrow4, ncol4, ...
                              nrow5, ncol5, ...
                              nvec );
disp(sprintf('total_mem_use=%g, mem_use5=%g', ...
              total_mem_use,    mem_use5 ));



Y = kronmult5(A1,A2,A3,A4,A5, X );
