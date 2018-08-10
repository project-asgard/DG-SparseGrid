function [x,error,iter,flag]=cg(x,b,max_it,tol)
%-------------------------------------------------------------
% cg.m solves the symmetric positive definite linear system Ax=b 
% using the Conjugate Gradient method with preconditioning.
%
% input   A        REAL symmetric positive definite matrix
%         x        REAL initial guess vector
%         b        REAL right hand side vector
%         M        REAL preconditioner matrix
%         max_it   INTEGER maximum number of iterations
%         tol      REAL error tolerance
%
% output  x        REAL solution vector
%         error    REAL error norm
%         iter     INTEGER number of iterations performed
%         flag     INTEGER: 0 = solution found to tolerance
%                           1 = no convergence given max_it
%---Modified by Lin Mu, 9/04/2015--------------------------

  flag = 0;                                 % initialization
  iter = 0;

  bnrm2 = norm( b );
  if  ( bnrm2 == 0.0 ), bnrm2 = 1.0; end
% full(Fc(x))
%  r = b - A*x;
%    r = b - Fc(x);

   r = b - ApplyA(x);

   
  error = norm( r ) / bnrm2;
  fprintf(' PCG residual(%d) = %g\n',iter,error)
  if ( error < tol ) return, end

  for iter = 1:max_it                       % begin iteration

     z  =r;%pc(r);% 
     rho = (r'*z);

     if ( iter > 1 ),                       % direction vector
        beta(iter) = rho / rho_1;
        p = z + beta(iter)*p;
     else
        p = z;
     end

%     q = A*p;
%       q = Fc(p);
      q = ApplyA(p);
      
     alpha(iter) = rho / (p'*q );
     x = x + alpha(iter) * p;                    % update approximation vector

     r = r - alpha(iter)*q;                      % compute residual
     error = norm( r ) / bnrm2;            % check convergence
     fprintf(' PCG residual(%d) = %g\n',iter,error)
     if ( error <= tol ), break, end 
       rho_1 = rho;

  end

%
%  compute eigenvalues and codition number
%
     d(1)= 1/alpha(1);
     for i = 2: iter
     d(i)=beta(i)/alpha(i-1)+1/alpha(i);
     end
     for i = 1: iter-1
     s(i)=-1*sqrt(beta(i+1))/alpha(i);
     end
%
     T = sparse(zeros(iter,iter));
     T(1,1)= d(1);
     for i=2:iter
     T(i,i) = d(i);
     end
     for i = 1:iter-1
     T(i,i+1)= s(i); T(i+1,i) = T(i,i+1);
     end
     lambda = eig(T);
%      lambda = eigs(T);
     lambdamax = max(lambda);
     lambdamin = min(lambda);
%      lambdamax = eigs(T,1,'LM');
%      lambdamin = eigs(T,1,'SM');
     condnumber = lambdamax/lambdamin;
     fprintf(' lambdamax = %f\n',lambdamax)
     fprintf(' lambdamin = %f\n',lambdamin)
     fprintf(' condnumber= %f\n',condnumber)

  if ( error > tol ) flag = 1; end         % no convergence
     fprintf(' PCG residual(%d) = %g\n',iter,error)

% END cg.m
  
 return;