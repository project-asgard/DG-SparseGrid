function tests = fk6d_test()

tests = functiontests(localfunctions);

end

function fk6d_vlasov4_fast_test(testCase)

disp('Testing with vlasov4 (fast implementation)');

quiet = 1;
lev = 3;
deg = 2;
TEND = 1;
slow = 0;
act_f = fk6d(Vlasov4,lev,deg,TEND,quiet,slow);

load('tests/vlasov4/solution.mat');

exp_f = fval;

verifyEqual(testCase,act_f,exp_f,'RelTol',1e-4)

end


function fk6d_vlasov4_slow_test(testCase)

disp('Testing with vlasov4 (slow reference implementation)');

quiet = 1;
lev = 3;
deg = 2;
TEND = 1;
slow = 1;
act_f = fk6d(Vlasov4,lev,deg,TEND,quiet,slow);

load('tests/vlasov4/solution.mat');

exp_f = fval;

verifyEqual(testCase,act_f,exp_f,'RelTol',1e-4)

end

