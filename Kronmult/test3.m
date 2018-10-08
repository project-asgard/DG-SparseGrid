
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

% ------------
[flops1,flops2,imethod] = flops_kron2( nrow1,ncol1, ...
                                       nrow2,ncol2 );
disp(sprintf('(A1)=(%d,%d),(A2)=(%d,%d),flops1=%g,flops2=%g', ...
              nrow1,ncol1,    nrow2,ncol2, flops1,flops2 ));

[flops1,flops2,imethod] = flops_kron2( nrow2,ncol2, ...
                                       nrow3,ncol3 );
disp(sprintf('(A2)=(%d,%d),(A3)=(%d,%d),flops1=%g,flops2=%g', ...
              nrow2,ncol2,    nrow3,ncol3, flops1,flops2 ));

% ------------

[flops1,flops2,imethod] = flops_kron3( nrow1,ncol1, ...
                                       nrow2,ncol2, ...
                                       nrow3,ncol3 );
disp(sprintf('(A1)=(%d,%d),(A2)=(%d,%d),(A3)=(%d,%d),flops1=%g,flops2=%g', ...
    nrow1,ncol1,  nrow2,ncol2,  nrow3,ncol3, flops1, flops2 ));


[flops1,flops2,imethod] = flops_kron3( nrow2,ncol2, ...
                                       nrow3,ncol3, ...
                                       nrow4,ncol4 );
disp(sprintf('(A2)=(%d,%d),(A3)=(%d,%d),(A4)=(%d,%d),flops1=%g,flops2=%g', ...
    nrow2,ncol2,  nrow3,ncol3,  nrow4,ncol4, flops1, flops2 ));


% ------------
[flops1,flops2,imethod] = flops_kron4( nrow1,ncol1, ...
                                       nrow2,ncol2, ...
                                       nrow3,ncol3, ...
                                       nrow4,ncol4 );
disp(sprintf('(A1)=(%d,%d),(A2)=(%d,%d),(A3)=(%d,%d),(A4)=(%d,%d),flops1=%g,flops2=%g', ...
    nrow1,ncol1,  nrow2,ncol2,  nrow3,ncol3, nrow4,ncol4, flops1, flops2 ));


[flops1,flops2,imethod] = flops_kron4( nrow2,ncol2, ...
                                       nrow3,ncol3, ...
                                       nrow4,ncol4, ...
                                       nrow5,ncol5 );
disp(sprintf('(A2)=(%d,%d),(A3)=(%d,%d),(A4)=(%d,%d),(A5)=(%d,%d),flops1=%g,flops2=%g', ...
    nrow2,ncol2,  nrow3,ncol3,  nrow4,ncol4, nrow5,ncol5, flops1, flops2 ));



% ------------
[flops1,flops2,imethod] = flops_kron5( nrow1,ncol1, ...
                                       nrow2,ncol2, ...
                                       nrow3,ncol3, ...
                                       nrow4,ncol4, ...
                                       nrow5,ncol5 );
disp(sprintf('(A1)=(%d,%d),(A2)=(%d,%d),(A3)=(%d,%d),(A4)=(%d,%d),(A5)=(%d,%d),flops1=%g,flops2=%g', ...
    nrow1,ncol1,  nrow2,ncol2,  nrow3,ncol3, nrow4,ncol4, nrow5,ncol5, flops1, flops2 ));


