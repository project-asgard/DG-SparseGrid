% simple script to test apply_FMWT
%
is_show_gflops = 0;

kdeg = 10;
Lev = 8;
n = kdeg * 2^Lev;
nvec = n;
X = rand(n,nvec);

tic();
FMWT =  OperatorTwoScale(kdeg,2^Lev);
time_fmwt = toc();
disp(sprintf('kdeg=%d,Lev=%d,n=%d, time for OperatorTwoScale is %g',...
              kdeg,   Lev,   n,  time_fmwt));
FMWT = full(FMWT);
% ------------------------------
% generate sparse matrix version
% ------------------------------
[ii,jj,aa] = find( FMWT );
FMWTs = sparse(ii,jj,aa,n,n);
clear ii jj aa;

imethod = 1;
tic;
Y_1 = apply_FMWT(kdeg,Lev, FMWT, X, imethod);
time_1 = toc();

imethod = 2;
tic;
Y_2 = apply_FMWT(kdeg,Lev, FMWT, X, imethod);
time_2 = toc();

imethod = 3;
tic;
Y_3 = apply_FMWT(kdeg,Lev, FMWT, X, imethod);
time_3 = toc();


err12 = norm(Y_1 - Y_2,1); 
relerr12 = err / max( norm(Y_1,1), norm(Y_2,1) );

err13 = norm(Y_1 - Y_3,1); 
relerr13 = err / max( norm(Y_1,1), norm(Y_3,1) );

disp(sprintf('kdeg=%d,Lev=%d,err12=%g,relerr12=%g, err13=%g, relerr13=%g',...
              kdeg,   Lev,   err12, relerr12,    err13,relerr13 ));

disp(sprintf('time_1=%g, time_2=%g, time_3=%g', ...
              time_1,    time_2,    time_3 ));

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
Y_s = FMWTs * X;
time_s = toc;
disp(sprintf('nnz(FMWTs)=%g, time_s=%g', ...
              nnz(FMWTs),    time_s ));

if (is_show_gflops),
  gflops = 2.0*nnz(FMWTs) * 1e-9 / time_s;
  disp(sprintf('sparse Gflops rate is %g Gflops/sec', ...
                gflops ));
end;

