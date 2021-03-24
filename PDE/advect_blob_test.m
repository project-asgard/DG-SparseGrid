function tests = advect_blob_test()
fh = localfunctions;
tests = functiontests(fh);
end

function advect_blob1_test(testCase)
addpath(genpath(pwd));
disp('Testing advect_blob1');
% setup PDE
args = {'lev',4,'deg',5,'dt',0.002,'quiet',true,'num_steps',5};
opts = OPTS(args);
pde = advect_blob1(opts);
% modify PDE - not needed here
% run PDE
[err,fval,fval_realspace,nodes,err_realspace] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = err / norm(fval);
verifyLessThan(testCase,rel_err,2e-4);
end

function advect_blob2_test(testCase)
addpath(genpath(pwd));
disp('Testing advect_blob2');
% setup PDE
args = {'lev',4,'deg',5,'dt',0.002,'quiet',true,'num_steps',5};
opts = OPTS(args);
pde = advect_blob2(opts);
% modify PDE - not needed here
% run PDE
[err,fval,fval_realspace,nodes,err_realspace] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = err / norm(fval);
verifyLessThan(testCase,rel_err,3e-4);
end