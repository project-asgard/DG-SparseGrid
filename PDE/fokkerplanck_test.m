function tests = fokkerplanck_test()
fh = localfunctions;
tests = functiontests(fh);
end

function fp1_C_pitch_test(testCase)
addpath(genpath(pwd));
disp('Testing the pitch dimension for C term within fokkerplanck1');
% setup PDE
args = {'lev',4,'deg',4,'dt',1e-3,'quiet',true,'num_steps',20,'timestep_method','matrix_exponential'};
opts = OPTS(args);
pde = fokkerplanck1_pitch_C(opts);
% modify PDE - not needed here
% run PDE
[err,fval,fval_realspace,nodes,err_realspace] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = err / norm(fval);
verifyLessThan(testCase,rel_err,1e-5);
end

function fp1_E_pitch_test(testCase)
addpath(genpath(pwd));
disp('Testing the pitch dimension for E term within fokkerplanck1');

% case 1 - flat initial condition

% setup PDE
args = {'lev',4,'deg',4,'dt',1e-1,'quiet',true,'num_steps',20,'timestep_method','matrix_exponential','case',1};
opts = OPTS(args);
pde = fokkerplanck1_pitch_E(opts);
% modify PDE - not needed here
% run PDE
[err,fval,fval_realspace,nodes,err_realspace] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = err / norm(fval);
verifyLessThan(testCase,rel_err,1e-4);

% case 2 - guassian initial condition

% setup PDE
args = {'lev',5,'deg',6,'dt',1e-1,'quiet',true,'num_steps',10,'timestep_method','matrix_exponential','case',2};
opts = OPTS(args);
pde = fokkerplanck1_pitch_E(opts);
% modify PDE - not needed here
% run PDE
[err,fval,fval_realspace,nodes,err_realspace] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = err / norm(fval);
verifyLessThan(testCase,rel_err,1e-4);
end

function fp1_R_pitch_test(testCase)
addpath(genpath(pwd));
disp('Testing the pitch dimension for R term within fokkerplanck1');

% setup PDE - case 1
args = {'lev',4,'deg',6,'dt',1e-1,'quiet',true,'num_steps',20,'timestep_method','matrix_exponential','case',1};
opts = OPTS(args);
pde = fokkerplanck1_pitch_R(opts);
% modify PDE - not needed here
% run PDE
[err,fval,fval_realspace,nodes,err_realspace] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = err / norm(fval);
verifyLessThan(testCase,rel_err,1e-5);

% setup PDE - case 2
args = {'lev',4,'deg',6,'dt',1e-1,'quiet',true,'num_steps',10,'timestep_method','matrix_exponential','case',2};
opts = OPTS(args);
pde = fokkerplanck1_pitch_R(opts);
% modify PDE - not needed here
% run PDE
[err,fval,fval_realspace,nodes,err_realspace] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = err / norm(fval);
verifyLessThan(testCase,rel_err,1e-4);

% setup PDE - case 3
args = {'lev',4,'deg',6,'dt',1e-1,'quiet',true,'num_steps',10,'timestep_method','matrix_exponential','case',3};
opts = OPTS(args);
pde = fokkerplanck1_pitch_R(opts);
% modify PDE - not needed here
% run PDE
[err,fval,fval_realspace,nodes,err_realspace] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = err / norm(fval);
verifyLessThan(testCase,rel_err,1e-4);

% setup PDE - case 4
args = {'lev',4,'deg',6,'dt',1e-1,'quiet',true,'num_steps',10,'timestep_method','matrix_exponential','case',4};
opts = OPTS(args);
pde = fokkerplanck1_pitch_R(opts);
% modify PDE - not needed here
% run PDE
[err,fval,fval_realspace,nodes,err_realspace] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = err / norm(fval);
verifyLessThan(testCase,rel_err,1e-4);

end
