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

function fp1_CE_pitch_test(testCase)
addpath(genpath(pwd));
disp('Testing the pitch dimension for C and E terms within fokkerplanck1');

% setup PDE - case 1
args = {'lev',4,'deg',6,'dt',1e-1,'quiet',true,'num_steps',20,'timestep_method','matrix_exponential','case',1};
opts = OPTS(args);
pde = fokkerplanck1_pitch_CE(opts);
% modify PDE - not needed here
% run PDE
[err,fval,fval_realspace,nodes,err_realspace] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = err / norm(fval);
verifyLessThan(testCase,rel_err,1e-5);

% setup PDE - case 2
args = {'lev',4,'deg',6,'dt',50e-1,'quiet',true,'num_steps',10,'timestep_method','matrix_exponential','case',2};
opts = OPTS(args);
pde = fokkerplanck1_pitch_CE(opts);
% modify PDE - not needed here
% run PDE
[err,fval,fval_realspace,nodes,err_realspace] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = err / norm(fval);
verifyLessThan(testCase,rel_err,1e-4);

% setup PDE - case 3
args = {'lev',4,'deg',6,'dt',50e-1,'quiet',true,'num_steps',1,'timestep_method','matrix_exponential','case',3};
opts = OPTS(args);
pde = fokkerplanck1_pitch_CE(opts);
% modify PDE - not needed here
% run PDE
[err,fval,fval_realspace,nodes,err_realspace,output] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = output.rel_err{end};
verifyLessThan(testCase,rel_err,1e-4);

% setup PDE - case 4
args = {'lev',4,'deg',6,'dt',50e-1,'quiet',true,'num_steps',1,'timestep_method','matrix_exponential','case',4};
opts = OPTS(args);
pde = fokkerplanck1_pitch_CE(opts);
% modify PDE - not needed here
% run PDE
[err,fval,fval_realspace,nodes,err_realspace,output] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = output.rel_err{end};
verifyLessThan(testCase,rel_err,1e-4);

end

function fp1_CER_pitch_test(testCase)
addpath(genpath(pwd));
disp('Testing the pitch dimension for C, E and R terms within fokkerplanck1');

% setup PDE - case 1
args = {'lev',4,'deg',7,'dt',100e-1,'quiet',true,'num_steps',1,'timestep_method','matrix_exponential','case',1,'normalize_by_mass',true};
opts = OPTS(args);
pde = fokkerplanck1_pitch_CER(opts);
% modify PDE - not needed here
% run PDE
[err,fval,fval_realspace,nodes,err_realspace,output] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = output.rel_err{end};
verifyLessThan(testCase,rel_err,1e-5);

end

function fp1_C_momentum_test(testCase)
addpath(genpath(pwd));
disp('Testing the momentum dimension for C terms within fokkerplanck1');

% setup PDE - case 1 (Maxwellian initial condition)
args = {'lev',4,'deg',4,'dt',1e2,'quiet',true,'num_steps',20,'timestep_method','matrix_exponential','case',1,'normalize_by_mass',false};
opts = OPTS(args);
pde = fokkerplanck1_momentum_C(opts);
% modify PDE - not needed here
% run PDE
[err,fval,fval_realspace,nodes,err_realspace,output] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = output.rel_err{end};
verifyLessThan(testCase,rel_err,2e-4);

% setup PDE - case 2 (step function initial condition)
args = {'lev',4,'deg',4,'dt',1e2,'quiet',true,'num_steps',20,'timestep_method','matrix_exponential','case',2,'normalize_by_mass',false};
opts = OPTS(args);
pde = fokkerplanck1_momentum_C(opts);
% modify PDE - not needed here
% run PDE
[err,fval,fval_realspace,nodes,err_realspace,output] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = output.rel_err{end};
verifyLessThan(testCase,rel_err,2e-4);

end

function fp1_C_momentum_LHS_test(testCase)
addpath(genpath(pwd));
disp('Testing the momentum dimension for C terms within fokkerplanck1_LHS');

% setup PDE - case 1 (Maxwellian initial condition)
args = {'lev',4,'deg',4,'dt',2,'quiet',true,'num_steps',50,'timestep_method','BE','case',1,'normalize_by_mass',false};
opts = OPTS(args);
pde = fokkerplanck1_momentum_C_LHS(opts);
% modify PDE - not needed here
% run PDE
[err,fval,fval_realspace,nodes,err_realspace,output] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = output.rel_err{end};
verifyLessThan(testCase,rel_err,2e-4);

% setup PDE - case 2 (step function initial condition)
args = {'lev',4,'deg',4,'dt',2,'quiet',true,'num_steps',50,'timestep_method','BE','case',2,'normalize_by_mass',false};
opts = OPTS(args);
pde = fokkerplanck1_momentum_C_LHS(opts);
% modify PDE - not needed here
% run PDE
[err,fval,fval_realspace,nodes,err_realspace,output] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = output.rel_err{end};
verifyLessThan(testCase,rel_err,3e-4);

end

function fp2_C_test(testCase)
addpath(genpath(pwd));
disp('Testing the momentum and pitch dimensions for term C within fokkerplanck2_C');

% setup PDE - case 1 (Step function initial condition)
args = {'lev',4,'deg',5,'dt',500,'quiet',true,'num_steps',1,'timestep_method','matrix_exponential','case',1,'normalize_by_mass',false};
opts = OPTS(args);
pde = fokkerplanck2_C(opts);
% modify PDE - not needed here
% run PDE
[err,fval,fval_realspace,nodes,err_realspace,output] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = output.rel_err{end};
verifyLessThan(testCase,rel_err,1e-4);

% setup PDE - case 2 (Maxwellian initial condition)
args = {'lev',4,'deg',5,'dt',500,'quiet',true,'num_steps',1,'timestep_method','matrix_exponential','case',2,'normalize_by_mass',false};
opts = OPTS(args);
pde = fokkerplanck2_C(opts);
% modify PDE - not needed here
% run PDE
[err,fval,fval_realspace,nodes,err_realspace,output] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = output.rel_err{end};
verifyLessThan(testCase,rel_err,1e-4);

end

function fp1_E_momentum_test(testCase)
addpath(genpath(pwd));
disp('Testing the momentum dimensions for term E within fokkerplanck1_momentum_E');

% setup PDE - case 1 (flat function initial condition)
args = {'lev',4,'deg',5,'dt',0.1,'quiet',true,'num_steps',1,'timestep_method','matrix_exponential','case',1,'normalize_by_mass',false};
opts = OPTS(args);
pde = fokkerplanck1_momentum_E(opts);
% modify PDE - not needed here
% run PDE
[err,fval,fval_realspace,nodes,err_realspace,output] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = output.rel_err{end};
verifyLessThan(testCase,rel_err,1e-4);

% setup PDE - case 2 (Maxwellian initial condition)
args = {'lev',4,'deg',5,'dt',0.1,'quiet',true,'num_steps',1,'timestep_method','matrix_exponential','case',2,'normalize_by_mass',false};
opts = OPTS(args);
pde = fokkerplanck1_momentum_E(opts);
% modify PDE - not needed here
% run PDE
[err,fval,fval_realspace,nodes,err_realspace,output] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = output.rel_err{end};
verifyLessThan(testCase,rel_err,1e-4);

end

function fp2_E_test(testCase)
addpath(genpath(pwd));
disp('Testing term E in 2D within fokkerplanck2_E');

% setup PDE - case 1 (flat function initial condition - solution does not change)
args = {'lev',4,'deg',4,'dt',0.01,'quiet',true,'num_steps',1,'timestep_method','matrix_exponential','case',1,'normalize_by_mass',false};
opts = OPTS(args);
pde = fokkerplanck2_E(opts);
% modify PDE - not needed here
% run PDE
[err,fval,fval_realspace,nodes,err_realspace,output] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = output.rel_err{end};
verifyLessThan(testCase,rel_err,1e-12);

% setup PDE - case 3 (manufactured solution)
args = {'lev',4,'deg',4,'dt',0.01,'quiet',false,'num_steps',1,'timestep_method','BE','case',3,'normalize_by_mass',false};
opts = OPTS(args);
pde = fokkerplanck2_E(opts);
% modify PDE - not needed here
% run PDE
[err,fval,fval_realspace,nodes,err_realspace,output] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = output.rel_err{end};
verifyLessThan(testCase,rel_err,1e-5);

end

function fp2_E_LHS_test(testCase)
addpath(genpath(pwd));
disp('Testing term E in 2D within fokkerplanck2_E');

% setup PDE - case 1 (flat function initial condition - solution does not change)
args = {'lev',4,'deg',4,'dt',0.01,'quiet',true,'num_steps',1,'timestep_method','matrix_exponential','case',1,'normalize_by_mass',false};
opts = OPTS(args);
pde = fokkerplanck2_E_LHS(opts);
% modify PDE - not needed here
% run PDE
[err,fval,fval_realspace,nodes,err_realspace,output] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = output.rel_err{end};
verifyLessThan(testCase,rel_err,1e-12);

% setup PDE - case 3 (manufactured solution)
args = {'lev',4,'deg',4,'dt',0.01,'quiet',false,'num_steps',1,'timestep_method','BE','case',3,'normalize_by_mass',false};
opts = OPTS(args);
pde = fokkerplanck2_E_LHS(opts);
% modify PDE - not needed here
% run PDE
[err,fval,fval_realspace,nodes,err_realspace,output] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = output.rel_err{end};
verifyLessThan(testCase,rel_err,1e-5);

end