function tests = projecti_test()
fh = localfunctions;
tests = functiontests(fh);
end

function pi_diff_test(testCase)
addpath(genpath(pwd));
disp('Testing the pitch dimension for C term within fokkerplanck1');
% setup PDE
args = {'lev',3,'deg',4,'dt',50,'quiet',true,'num_steps',10,'timestep_method','matrix_exponential','case',1};
opts = OPTS(args);
pde = projecti_diff1(opts);
% modify PDE - not needed here
% run PDE
[err,fval,fval_realspace,nodes,err_realspace] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = err / norm(fval);
verifyLessThan(testCase,rel_err,1e-5);

end