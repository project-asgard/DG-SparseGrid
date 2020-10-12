function tests = moment_integral_test()
tests = functiontests(localfunctions);
end

function test_moment_homo_1D(testCase)
%1D Homogeneous case
addpath(genpath(pwd));
[err, act_f, act_frs] = asgard(advection1, 'quiet', true, 'num_steps', 5, 'deg', 3, 'lev', 3);
test_func = @(x,p,t) x.*0 + 1;
dim_x.domainMin = 0;
dim_x.domainMax = pi;
dim_x = check_dimension(1,dim_x);
dimensions = {dim_x};
test_moment_zero = moment_integral(3, 3, act_frs, test_func, dimensions);
verifyLessThan(testCase, test_moment_zero, 1e-6);
end

function test_moment_inhomo_1D(testCase)
%1D Inhomogeneous Case
addpath(genpath(pwd));
[err, act_f, act_frs] = asgard(advection1, 'quiet', true, 'num_steps', 5, 'deg', 3, 'lev', 3);
dim_x.domainMin = 0;
dim_x.domainMax = pi;
dim_x = check_dimension(1,dim_x);
dimensions = {dim_x};
test_func_nonzero = @(x,p,t) cos(x);
test_moment_nonzero = moment_integral(3, 3, act_frs, test_func_nonzero, dimensions);
diff = abs(test_moment_nonzero - pi/2);
verifyLessThan(testCase, diff, 1e-6);
end

function test_moment_2D(testCase)
[err, act_f, act_frs] = asgard(continuity2, 'timestep_method','CN','quiet', true, 'num_steps', 5, 'deg', 3, 'lev',3);
dim_x.domainMin = -1;
dim_x.domainMax = +1;
dim_y.domainMin = -2;
dim_y.domainMax = +2;
test_func = @(p,z,t) p.*0 + 1;
dim_x = check_dimension(2,dim_x);
dim_y = check_dimension(2,dim_y);
dimensions = {dim_x, dim_y};
test_moment_2D = moment_integral(3, 3, act_frs, test_func, dimensions);
diff = abs(test_moment_2D);
verifyLessThan(testCase, diff, 1e-6);
end