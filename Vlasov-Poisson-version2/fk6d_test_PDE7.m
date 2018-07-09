function tests = fk6d_test_PDE7()
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

  if (isOctave),
    % -----------------------------
    % perform octave test framework
    % -----------------------------
    reporting_level = 'quiet';
    reporting_level = 'verbose';
    reporting_level = 'normal';
  
%     [n,nmax,nxfail,nskip] = test( 'fk6d_vlasov4_check', reporting_level);
    [n,nmax,nxfail,nskip] = test( 'fk6d_vlasov7_check', reporting_level);
    disp(sprintf('%d tests were executed', nmax));
    disp(sprintf('%d tests passed ',n));
    disp(sprintf('%d tests failed ', nxfail));
    disp(sprintf('%d tests were skipped', nskip));
  
  
  else
    fh = localfunctions;
    tests = functiontests(fh);
  end;


end

function fk6d_vlasov7_compression0_test(testCase)

addpath(genpath(pwd));

disp('Testing with vlasov7 (compression==0)');

quiet = 1;
lev = 3;
deg = 2;
TEND = 1;
compression = 0;
act_f = fk6d(Vlasov7,lev,deg,TEND,quiet,compression);

load('tests/vlasov7/solution.mat');

exp_f = fval;

verifyEqual(testCase,act_f,exp_f,'RelTol',1e-4);

end


function fk6d_vlasov7_compression1_test(testCase)

addpath(genpath(pwd));

disp('Testing with vlasov7 (compression==1)');

quiet = 1;
lev = 3;
deg = 2;
TEND = 1;
compression = 1;
act_f = fk6d(Vlasov7,lev,deg,TEND,quiet,compression);

load('tests/vlasov7/solution.mat');

exp_f = fval;

verifyEqual(testCase,act_f,exp_f,'RelTol',1e-4);

end

function fk6d_vlasov7_compression2_test(testCase)

addpath(genpath(pwd));

disp('Testing with vlasov7 (compression==2)');

quiet = 1;
lev = 3;
deg = 2;
TEND = 1;
compression = 2;
act_f = fk6d(Vlasov7,lev,deg,TEND,quiet,compression);

load('tests/vlasov7/solution.mat');

exp_f = fval;

verifyEqual(testCase,act_f,exp_f,'RelTol',1e-4);

end

function fk6d_vlasov7_compression3_test(testCase)

addpath(genpath(pwd));

disp('Testing with vlasov7 (compression==3)');

quiet = 1;
lev = 3;
deg = 2;
TEND = 1;
compression = 3;
act_f = fk6d(Vlasov7,lev,deg,TEND,quiet,compression);

load('tests/vlasov7/solution.mat');

exp_f = fval;

verifyEqual(testCase,act_f,exp_f,'RelTol',1e-4);

end

function fk6d_vlasov7_compression4_test(testCase)

addpath(genpath(pwd));

disp('Testing with vlasov7 (compression==4)');

quiet = 1;
lev = 3;
deg = 2;
TEND = 1;
compression = 4;
act_f = fk6d(Vlasov7,lev,deg,TEND,quiet,compression);

load('tests/vlasov7/solution.mat');

exp_f = fval;

verifyEqual(testCase,act_f,exp_f,'RelTol',1e-4);

end



