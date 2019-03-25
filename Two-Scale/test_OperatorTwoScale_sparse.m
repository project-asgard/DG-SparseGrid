% simple test for OperatorTwoScale_sparse
%
k = 8;
L = 4;
n = k*2^L;
x = linspace(0,1,n);

FMWTs = OperatorTwoScale_sparse( k, 2^L );
FMWT  = OperatorTwoScale(k, 2^L);

% -----------------------------------
% FMWT should be an orthogonal matrix
% -----------------------------------
errs = norm( FMWTs' * FMWTs - eye(n,n), 1 );
err  = norm( FMWT' * FMWT - eye(n,n), 1 );
disp(sprintf('orthonormal err=%g, sparse err=%g', ...
	                  err,    errs ));



figure(1);
subplot(2,1,1); spy( FMWT ); title('FMWT');
subplot(2,1,2); spy( FMWTs ); title('FMWTs');

for j=0:(k-1),
    Y  = FMWT  * reshape(x.^j,  n,1);
    Ys = FMWTs * reshape(x.^j,  n,1);
    kp1 = k+1;
    disp(sprintf('j=%d,norm(Y(kp1:n))=%g, norm(Ys(kp1:n))=%g', ...
		  j,   norm(Y(kp1:n)),    norm(Ys(kp1:n))  ));

    do_plot_Ys = 0;
    if (do_plot_Ys),
     figure();
     plot( kp1:n, Y(kp1:n),'rx', kp1:n, Ys(kp1:n), 'bo'); 
     title(sprintf('j=%d,norm(Y(kp1:n))=%g,norm(Ys(kp1:n))=%g',...
		    j,   norm(Y(kp1:n)),   norm(Ys(kp1:n))  ));
    end;
end;


figure();
subplot(2,1,1); plot( FMWT(1:k,:)'); title('FMWT');
subplot(2,1,2); plot( FMWTs(1:k,:)'); title('FMWTs');


do_plot_mesh = 0;
if (do_plot_mesh),
figure();
subplot(2,1,1); mesh( full(FMWT)); title('FMWT');
subplot(2,1,2); mesh( full(FMWTs)); title('FMWTs');
end;
