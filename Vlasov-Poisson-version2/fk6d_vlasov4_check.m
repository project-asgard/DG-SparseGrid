function isok = fk6d_vlasov4_check(testCase)
% isok = fk6d_vlasov4_check(testCase)
%
  addpath(genpath(pwd));
  disp(sprintf('Testing with valsov4, testCase=%d',testCase));
  
  switch(testCase)
   case 1 
      quiet = 1; lev = 3; deg = 2; TEND = 1; compression = 0;
   case 2 
      quiet = 1; lev = 3; deg = 2; TEND = 1; compression = 1;
   case 3 
      quiet = 1; lev = 3; deg = 2; TEND = 1; compression = 2;
   case 4 
      quiet = 1; lev = 3; deg = 2; TEND = 1; compression = 3;
   case 5 
      quiet = 1; lev = 3; deg = 2; TEND = 1; compression = 4;
   otherwise 
      error(sprintf('fk6d_vlasov4_check: invalid testCase=%d',testCase));
  end;
  
  act_f = fk6d(Vlasov4,lev,deg,TEND,quiet,compression);

  load('tests/vlasov4/solution.mat');
  exp_f = fval;
  
  isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
  
  tol = 1e-4;
  isok = all(   abs(act_f(:)-exp_f(:)) <= tol * abs( exp_f(:) ) );
  if (~isOctave),
    verifyEqual(testCase,act_f,exp_f,'RelTol',tol);
  end;

end 

%!test assert(fk6d_vlasov4_check(1))
%!test assert(fk6d_vlasov4_check(2))
%!test assert(fk6d_vlasov4_check(3))
%!test assert(fk6d_vlasov4_check(4))
%!test assert(fk6d_vlasov4_check(5))
