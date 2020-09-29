function tests = moment_integral_test()
tests = functiontests(localfunctions);
end

function testall(testCase)

addpath(genpath(pwd));
[err, act_f, act_frs] = asgard(advection1, 'quiet', true, 'num_steps', 5, 'deg', 3, 'lev', 3);
test_func = @(x,p,t) x.*0 + 1;
dim_x.domainMin = 0;
dim_x.domainMax = pi;
dimensions = {dim_x};
test_moment_zero = moment_integral(3, 3, act_frs, test_func, dimensions);
verifyLessThan(testCase, test_moment_zero, 1e-6);
test_func_nonzero = @(x,p,t) cos(x);
test_moment_nonzero = moment_integral(3, 3, act_frs, test_func_nonzero, dimensions);
diff = abs(test_moment_nonzero - pi/2);
verifyLessThan(testCase, diff, 1e-6);

end