function [flops,isplit,imethod] = kron_minflops( rc )
% -----------------------------------------------
% return min number of flops to evaluate kron product
% kron(A1, ..., Ak) where Aj is rc(1,j) by rc(2,j) 
% -----------------------------------------------
global kron_minflops_;
global kron_split_;
global kron_method_;
idebug = 1;

c = who('kron_minflops_');
is_first = (length(c) == 0);
if (is_first),
  rc0 = [1,1];
  key = sprintf('i%d_',rc0(:));
  kron_minflops_.(key) = 1;
  kron_split_.(key) = 1;
  kron_method_.(key) = 1;
end;

n = size(rc,2);
key = sprintf('i%d_',rc(:));
if (n == 1),
   nrow = rc(1,1);
   ncol = rc(2,1);
   flops = 2.0*nrow*ncol;
   isplit = 1;
   kron_minflops_.(key) = flops;
   kron_split_.(key) = 1;
   kron_method_.(key) = 1;

   flops = kron_minflops_.(key);
   isplit = kron_split_.(key);
   imethod = kron_method_.(key);

   return;
end;

if (isfield(kron_minflops_,key)),
   flops = kron_minflops_.(key);
   isplit = kron_split_.(key);
   imethod = kron_method_.(key);
   return;
end;

% ------------------------------------------
% new entry, need to explore different cases
% ------------------------------------------
for k=1:(n-1),
  % -------------------------- 
  % Y = kron(A1,..An) * X
  % Y = kron(Akp1..An) * X * kron(A1..Ak))^t
  % 
  %  m1 = prod( rc(1, 1:k) );
  %  m2 = prod( rc(1, kp1:n) );
  %  n2 = prod( rc(2, kp1:n) ); 
  %  n1 = prod( rc(2, 1:k) );
  %
  % ---------
  % method 1:
  % ---------
  %  X = reshape( X, n2,n1);
  %  T1 = kron(Akp1..An)*X
  %  T1 is m2 by n1 
  %  T1t = reshape( transpose(T1), n1, m2)
  %  Y = (kron(A1..Ak) T1t)^t
  %  
  %
  % ---------
  % method 2:
  % ---------
  % X = reshape( X, n2, n1 )
  % T2 = kron(A1...Ak)* X^t
  % T2 is  m1 by n2
  % T2t = reshape( transpose(T2), n2, m1)
  % Y = kron(Akp1..An) * T2t
  % -------------------------- 
  kp1 = k + 1;
  n1 = prod( rc(2,1:k) );
  n2 = prod( rc(2,kp1:n));
  m1 = prod( rc(1,1:k));
  m2 = prod( rc(1,kp1:n));

  % --------
  % method 1
  % --------
  %  --------------------------------------
  %  X = reshape( X, n2,n1);
  %  T1 = kron(Akp1..An)*X
  %  T1 is m2 by n1 
  %  T1t = reshape( transpose(T1), n1, m2)
  %  Y = (kron(A1..Ak) T1t)^t
  %  --------------------------------------
  [flops1,jsplit1,jmethod1] = kron_minflops( rc(:,kp1:n) );
  flops_t1 =   flops1 *  n1;
  [flops2,jsplit2,jmethod2] = kron_minflops( rc(:,1:k) );
  flops_method_1 = flops2 * m2 + flops_t1;
   
  % --------
  % method 2
  % --------
  %  --------------------------------------
  % X = reshape( X, n1, n2 )
  % T2 = kron(A1...Ak)* X^t
  % T2 is  m1 by n2
  % T2t = reshape( transpose(T2), n2, m1)
  % Y = kron(Akp1..An) * T2t
  %  --------------------------------------
  [flops1,jsplit1] = kron_minflops( rc(:,1:k) );
  flops_t2 = flops1 * n2;
  [flops2,jsplit2] = kron_minflops( rc(:,kp1:n) );
  flops_method_2 = flops2 * m1 + flops_t2;
  if (idebug >= 1),
    disp(sprintf('kron_minflops:n=%d,k=%d,flops_method_1=%g, flops_method_2=%g',...
         n, k, flops_method_1, flops_method_2 ));
  end;
  if (flops_method_1 <= flops_method_2),
     total_flops(k) = flops_method_1;
     kron_method(k) = 1;
  else
     total_flops(k) = flops_method_2;
     kron_method(k) = 2;
  end;
end;
[flops, isplit] = min( total_flops );
imethod = kron_method(isplit);

kron_minflops_.(key) = flops;
kron_split_.(key) = k;
kron_method_.(key) = imethod;


