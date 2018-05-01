function tests = fk6d_test()

tests = functiontests(localfunctions);

end

function fk6d_vlasov4_compression0_test(testCase)

addpath(genpath(pwd));

disp('Testing with vlasov4 (compression==0)');

quiet = 1;
lev = 3;
deg = 2;
TEND = 1;
compression = 0;
act_f = fk6d(Vlasov4,lev,deg,TEND,quiet,compression);

load('tests/vlasov4/solution.mat');

exp_f = fval;

verifyEqual(testCase,act_f,exp_f,'RelTol',1e-4);

end


function fk6d_vlasov4_compression1_test(testCase)

addpath(genpath(pwd));

disp('Testing with vlasov4 (compression==1)');

quiet = 1;
lev = 3;
deg = 2;
TEND = 1;
compression = 1;
act_f = fk6d(Vlasov4,lev,deg,TEND,quiet,compression);

load('tests/vlasov4/solution.mat');

exp_f = fval;

verifyEqual(testCase,act_f,exp_f,'RelTol',1e-4);

end

function fk6d_vlasov4_compression2_test(testCase)

addpath(genpath(pwd));

disp('Testing with vlasov4 (compression==2)');

quiet = 1;
lev = 3;
deg = 2;
TEND = 1;
compression = 2;
act_f = fk6d(Vlasov4,lev,deg,TEND,quiet,compression);

load('tests/vlasov4/solution.mat');

exp_f = fval;

verifyEqual(testCase,act_f,exp_f,'RelTol',1e-4);

end

function fk6d_vlasov4_compression3_test(testCase)

addpath(genpath(pwd));

disp('Testing with vlasov4 (compression==3)');

quiet = 1;
lev = 3;
deg = 2;
TEND = 1;
compression = 3;
act_f = fk6d(Vlasov4,lev,deg,TEND,quiet,compression);

load('tests/vlasov4/solution.mat');

exp_f = fval;

verifyEqual(testCase,act_f,exp_f,'RelTol',1e-4);

end

function fk6d_vlasov4_compression4_test(testCase)

addpath(genpath(pwd));

disp('Testing with vlasov4 (compression==4)');

quiet = 1;
lev = 3;
deg = 2;
TEND = 1;
compression = 4;
act_f = fk6d(Vlasov4,lev,deg,TEND,quiet,compression);

load('tests/vlasov4/solution.mat');

exp_f = fval;

verifyEqual(testCase,act_f,exp_f,'RelTol',1e-4);

end

