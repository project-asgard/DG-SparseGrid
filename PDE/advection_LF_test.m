function tests = advection1_LF_test()
fh = localfunctions;
tests = functiontests(fh);
end

function advection1_LF_left_test(testCase)
addpath(genpath(pwd));
disp('Testing advection1_LF the left(-1) flux');
% setup PDE
args = {'lev',4,'deg',4,'dt',1e-2,'quiet',true,'num_steps',20,'timestep_method','matrix_exponential'};
opts = OPTS(args);
pde = advection1_LF(opts);
% modify PDE - not needed here
% run PDE
[err,fval,fval_realspace,nodes,err_realspace] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = err / norm(fval);
verifyLessThan(testCase,rel_err,1e-3);
end

function advection2_LF_left_test(testCase)
addpath(genpath(pwd));
disp('Testing advection2_LF the left(-1) flux');
% setup PDE
args = {'lev',4,'deg',4,'dt',1e-2,'quiet',true,'num_steps',15,'timestep_method','BE'};
opts = OPTS(args);
pde = advection2_LF(opts);
% modify PDE - not needed here
% run PDE
[err,fval,fval_realspace,nodes,err_realspace] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = err / norm(fval);
verifyLessThan(testCase,rel_err,8e-3);
end

