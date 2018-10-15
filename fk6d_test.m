function tests = fk6d_test()
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

  if (isOctave),
    % -----------------------------
    % perform octave test framework
    % -----------------------------
    reporting_level = 'quiet';
    reporting_level = 'verbose';
    reporting_level = 'normal';
  
    [n,nmax,nxfail,nskip] = test( 'fk6d_vlasov4_check', reporting_level);
    disp(sprintf('%d tests were executed', nmax));
    disp(sprintf('%d tests passed ',n));
    disp(sprintf('%d tests failed ', nxfail));
    disp(sprintf('%d tests were skipped', nskip));
  
  
  else
    fh = localfunctions;
    tests = functiontests(fh);
  end;


end

function fk6d_vlasov4_compression0_test(testCase)

addpath(genpath(pwd));

disp('Testing with vlasov4 (compression==0)');

quiet = 1;
lev = 3;
deg = 2;
TEND = 1;
compression = 0;
[err,act_f,act_frs] = fk6d(Vlasov4,lev,deg,TEND,quiet,compression);

load('tests/vlasov4/solution.mat');

exp_f = fval;
exp_frs = fval_realspace;

%verifyEqual(testCase,act_f,exp_f,'RelTol',1e-4);
verifyEqual(testCase,act_frs,exp_frs,'RelTol',1e-4);

end


function fk6d_vlasov4_compression1_test(testCase)

addpath(genpath(pwd));

disp('Testing with vlasov4 (compression==1)');

quiet = 1;
lev = 3;
deg = 2;
TEND = 1;
compression = 1;
[err,act_f,act_frs] = fk6d(Vlasov4,lev,deg,TEND,quiet,compression);

load('tests/vlasov4/solution.mat');

exp_f = fval;
exp_frs = fval_realspace;

%verifyEqual(testCase,act_f,exp_f,'RelTol',1e-4);
verifyEqual(testCase,act_frs,exp_frs,'RelTol',1e-4);

end

function fk6d_vlasov4_compression2_test(testCase)

addpath(genpath(pwd));

disp('Testing with vlasov4 (compression==2)');

quiet = 1;
lev = 3;
deg = 2;
TEND = 1;
compression = 2;
[err,act_f,act_frs] = fk6d(Vlasov4,lev,deg,TEND,quiet,compression);

load('tests/vlasov4/solution.mat');

exp_f = fval;
exp_frs = fval_realspace;

%verifyEqual(testCase,act_f,exp_f,'RelTol',1e-4);
verifyEqual(testCase,act_frs,exp_frs,'RelTol',1e-4);

end

function fk6d_vlasov4_compression3_test(testCase)

addpath(genpath(pwd));

disp('Testing with vlasov4 (compression==3)');

quiet = 1;
lev = 3;
deg = 2;
TEND = 1;
compression = 3;
[err,act_f,act_frs] = fk6d(Vlasov4,lev,deg,TEND,quiet,compression);

load('tests/vlasov4/solution.mat');

exp_f = fval;
exp_frs = fval_realspace;

%verifyEqual(testCase,act_f,exp_f,'RelTol',1e-4);
verifyEqual(testCase,act_frs,exp_frs,'RelTol',1e-4);

end

function fk6d_vlasov4_compression4_test(testCase)

addpath(genpath(pwd));

disp('Testing with vlasov4 (compression==4)');

quiet = 1;
lev = 3;
deg = 2;
TEND = 1;
compression = 4;
[err,act_f,act_frs] = fk6d(Vlasov4,lev,deg,TEND,quiet,compression);

load('tests/vlasov4/solution.mat');

exp_f = fval;
exp_frs = fval_realspace;

%verifyEqual(testCase,act_f,exp_f,'RelTol',1e-4);
verifyEqual(testCase,act_frs,exp_frs,'RelTol',1e-4);

end

function fk6d_vlasov7_compression4_lev3_test(testCase)

addpath(genpath(pwd));

disp('Testing with vlasov7 (compression==4, Lev=3)');

quiet = 1;
lev = 3;
deg = 2;
TEND = 0.05;
compression = 4;
[err,act_f,act_frs] = fk6d(Vlasov7,lev,deg,TEND,quiet,compression);

tol = 1e-2;
assert(err(1) < tol);

end

function fk6d_vlasov7_compression4_lev4_test(testCase)

addpath(genpath(pwd));

disp('Testing with vlasov7 (compression==4, Lev=4)');

quiet = 1;
lev = 4;
deg = 2;
TEND = 0.05;
compression = 4;
[err,act_f,act_frs] = fk6d(Vlasov7,lev,deg,TEND,quiet,compression);

tol = 1e-2;
assert(err(1) < tol);

end

function fk6d_vlasov7_implicit_BE_test(testCase)

addpath(genpath(pwd));

disp('Testing with vlasov7 (compression==1, Lev=3, implicit==1)');

quiet = 1;
lev = 3;
deg = 2;
TEND = 0.05;
compression = 4;
implicit = 1;
[err,act_f,act_frs] = fk6d(Vlasov7,lev,deg,TEND,quiet,compression,implicit);

tol = 1e-2;
assert(err(1) < tol);

end

%!test fk6d_vlasov4_compression0_test(1)
%!test fk6d_vlasov4_compression1_test(2)
%!test fk6d_vlasov4_compression2_test(3)
%!test fk6d_vlasov4_compression3_test(4)
%!test fk6d_vlasov4_compression4_test(5)

