% simple test for OperatorTwoScale_wavelet2
%
maxDeg = 8;
nLev = 9;

k = maxDeg;
L = nLev;
n = k*2^L;

disp(sprintf('maxDeg=%d, nLev=%d, n = %d', ...
	      maxDeg,    nLev,    n ));

tic();
FMWT2 = OperatorTwoScale_wavelet2( k, 2^L );
time_wavelet2 = toc();

tic();
FMWT  = OperatorTwoScale_wavelet(k, 2^L);
time_wavelet = toc();



disp(sprintf('time for OperatorTwoScale_wavelet2=%g ', ...
	      time_wavelet2 ));
disp(sprintf('time for OperatorTwoScale_wavelet=%g ', ...
	      time_wavelet ));

disp(sprintf('norm 1 difference is %g ', ...
	      norm( FMWT2 - FMWT,1) ));

max_diff = max( abs( FMWT2(:) - FMWT(:) ) );
max_diff = full(max_diff);
disp(sprintf('max difference is %g ', ...
	      max_diff ));

disp(sprintf('norm(eye - trans(FMWT2)*FMWT2,1)=%g', ...
	      norm(speye(n,n) - FMWT2'*FMWT2,1) ));
disp(sprintf('norm(eye - trans(FMWT)*FMWT,1)=%g', ...
	      norm(speye(n,n) - FMWT'*FMWT, 1)  ));
