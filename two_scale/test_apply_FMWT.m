% simple script to test apply_FMWT
%
is_show_gflops = 0;

nerrors = 0;
tol = 1e-9;

kdeg = 4;
Lev = 8;
n = kdeg * 2^Lev;
nvec = n;
X = rand(n,nvec);

global OperatorTwoScale_method;

for icase=1:3,
  OperatorTwoScale_method = cmerge( icase==1,'wavelet', ...
				    cmerge(icase==2,'wavelet2', 'nonwavelet'));

  tic();
  FMWT = OperatorTwoScale(kdeg,Lev);
  time_fmwt = toc();
  disp('==========');
  if (strcmp(OperatorTwoScale_method,'wavelet')),
    disp(sprintf('kdeg=%d,Lev=%d,n=%d, time for OperatorTwoScale_wavelet is %g',...
                  kdeg,   Lev,   n,  time_fmwt));
  elseif (strcmp(OperatorTwoScale_method,'nonwavelet')),
    disp(sprintf('kdeg=%d,Lev=%d,n=%d, time for OperatorTwoScale_nonwavelet is %g',...
                  kdeg,   Lev,   n,  time_fmwt));
  
  elseif (strcmp(OperatorTwoScale_method,'wavelet2')),
    disp(sprintf('kdeg=%d,Lev=%d,n=%d, time for OperatorTwoScale_wavelet2 is %g',...
                  kdeg,   Lev,   n,  time_fmwt));
  end;
  disp('==========');



FMWT = full(FMWT);
% ------------------------------
% generate sparse matrix version
% ------------------------------

FMWTs = sparsify_matrix( FMWT );


for isLeft=1:-1:0,
for isTrans=0:1,

   LT_in = strcat( ...
             cmerge( isLeft,'L','R'), ...
             cmerge( isTrans,'T','N') );
   
   disp('------------');
   disp(sprintf('LT_in=%s', LT_in));
   disp('------------');


  imethod = 1;
  tic;
  Y_1 = apply_FMWT(kdeg,Lev, FMWT, X, imethod,LT_in);
  time_1 = toc();
  
  imethod = 2;
  tic;
  Y_2 = apply_FMWT(kdeg,Lev, FMWT, X, imethod,LT_in);
  time_2 = toc();
  
  imethod = 3;
  tic;
  Y_3 = apply_FMWT(kdeg,Lev, FMWT, X, imethod,LT_in);
  time_3 = toc();
  
  
  err12 = norm(Y_1 - Y_2,1); 
  relerr12 = err12 / max( norm(Y_1,1), norm(Y_2,1) );
  
  err13 = norm(Y_1 - Y_3,1); 
  relerr13 = err13 / max( norm(Y_1,1), norm(Y_3,1) );
  
  disp(sprintf('kdeg=%d,Lev=%d,err12=%g,relerr12=%g, err13=%g, relerr13=%g',...
                kdeg,   Lev,   err12, relerr12,    err13,relerr13 ));
  
  disp(sprintf('time_1=%g, time_2=%g, time_3=%g', ...
                time_1,    time_2,    time_3 ));
  
  if (relerr12 >= tol),
     nerrors = nerrors + 1;
  end;
  if (relerr13 >= tol),
     nerrors = nerrors + 1;
  end;

  if (is_show_gflops),
    gflops = 2.0*n*n * size(X,2) * 1e-9/time_1;
    disp(sprintf('Dense Gflops is %g Gflops/sec ', ...
                  gflops ));
    
    
    gflops = 2.0*nnz(FMWTs) * size(X,2) * 1e-9/time_2;
    disp(sprintf('Structure Dense Gflops is %g Gflops/sec ', ...
                  gflops ));
  
    gflops = 2.0*nnz(FMWTs) * size(X,2) * 1e-9/time_3;
    disp(sprintf('Structure Identical Dense Gflops is %g Gflops/sec ', ...
                  gflops ));
  end;
  
  % --------------------
  % test sparse version
  % --------------------
  
  tic;
  if (isLeft),
     Y_s = cmerge(isTrans, FMWTs',FMWTs) * X;
  else
     Y_s = X * cmerge(isTrans, FMWTs', FMWTs);
  end;
  time_s = toc;
  disp(sprintf('nnz(FMWTs)=%g, time_s=%g', ...
                nnz(FMWTs),    time_s ));
  
  if (is_show_gflops),
    gflops = 2.0*nnz(FMWTs) * 1e-9 / time_s;
    disp(sprintf('sparse Gflops rate is %g Gflops/sec', ...
                  gflops ));
  end;
  
  err1s = norm(Y_1 - Y_s,1); 
  relerr1s = err1s / max( norm(Y_1,1), norm(Y_s,1) );

  if (relerr1s >= tol),
    nerrors = nerrors + 1;
  end;
  
  disp(sprintf('kdeg=%d,Lev=%d,err1s=%g,relerr1s=%g ', ...
                kdeg,   Lev,   err1s,   relerr1s ));

end; % for isTrans
end; % for isLeft



end; % icase

if (nerrors == 0),
   disp('ALL OK');
else
  disp(sprintf('there are %d errors ', nerrors ));
end;

