function tests = asgard_test()

fh = localfunctions;
tests = functiontests(fh);

end

function asgard_advection1_explicit_test(testCase)

addpath(genpath(pwd));

disp('Testing advection1 (RK3)');

[err, act_f, act_frs] = asgard(advection1, 'quiet', true, 'num_steps', 5);

verifyLessThan(testCase, err, 9e-4);

end

function asgard_advection1reverse_explicit_test(testCase)

addpath(genpath(pwd));

disp('Testing advection1_reverse_flow (RK3)');

[err, act_f, act_frs] = asgard(advection1_reverse_flow, 'quiet', true, 'num_steps', 5);

verifyLessThan(testCase, err, 9e-4);

end

function asgard_advection1_implicit_test(testCase)

addpath(genpath(pwd));

disp('Testing advection1 (CN)');

[err, act_f, act_frs] = asgard(advection1, 'quiet',true,'timestep_method','CN','lev',4,'deg',3,'num_steps',5);

verifyLessThan(testCase, err, 3e-6);

end

function asgard_advection1reverse_implicit_test(testCase)

addpath(genpath(pwd));

disp('Testing advection1_reverse_flow (CN)');

[err, act_f, act_frs] = asgard(advection1_reverse_flow, ...
    'quiet',true,'timestep_method', 'CN','lev',4,'deg',3,'num_steps',5);

verifyLessThan(testCase, err, 3e-6);

end

function asgard_diffusion1_explicit_test(testCase)

addpath(genpath(pwd));

disp('Testing diffusion1 (RK3)');

[err,act_f,act_frs] = asgard(diffusion1,'lev',3,'quiet',true,'deg',3);

verifyLessThan(testCase,err,2.8e-5);

end

function asgard_diffusion1_implicit_test(testCase)

addpath(genpath(pwd));

disp('Testing diffusion1 (CN)');

[err,act_f,act_frs] = asgard(diffusion1,'lev',3,'quiet',true,'deg',3,'timestep_method', 'CN');

verifyLessThan(testCase,err,2.8e-5);

end


function asgard_continuity1_ex_test(testCase)

addpath(genpath(pwd));

disp('Testing continuity1 (RK3)');

[err,act_f,act_frs] = asgard(continuity1,'quiet',true,'num_steps',2);

verifyLessThan(testCase,err,1e-4);

end


function asgard_continuity1_adapt_test(testCase)

addpath(genpath(pwd));

disp('Testing continuity1 (CN/adapt)');

[err,act_f,act_frs] = asgard(continuity1,'quiet',true,...
    'timestep_method', 'CN','num_steps',2,'adapt',true);

verifyLessThan(testCase,err,1e-4);

end

function asgard_continuity1_oldhash_test(testCase)

addpath(genpath(pwd));

disp('Testing continuity1 (RK3/old_hash)');

[err,act_f,act_frs] = asgard(continuity1,'quiet',true,'num_steps',2,'use_oldhash',true);

verifyLessThan(testCase,err,1e-4);

end

function asgard_continuity2_expliciti_test(testCase)

addpath(genpath(pwd));

disp('Testing continuity2 (RK3)');

[err,act_f,act_frs] = asgard(continuity2,'quiet',true,'lev',4,'deg',3,'num_steps',2);

verifyLessThan(testCase,err,1e-4);

end

function asgard_continuity2_implicit_test(testCase)

addpath(genpath(pwd));

disp('Testing continuity2 (CN)');

[err,act_f,act_frs] = asgard(continuity2,'quiet',true,...
    'timestep_method', 'CN','lev',4,'deg',3,'num_steps',2);

verifyLessThan(testCase,err,1e-4);

end

function asgard_continuity2_adapt_test(testCase)

addpath(genpath(pwd));

disp('Testing continuity2 (CN/adapt)');

[err,act_f,act_frs] = asgard(continuity2,'quiet',true,...
    'timestep_method', 'CN','lev',4,'deg',3,'num_steps',2,'adapt',true);

verifyLessThan(testCase,err,1e-4);

end

function asgard_fokkerplanck4p1b_explicit_test(testCase)

addpath(genpath(pwd));

disp('Testing fokkerplanck4p1b (RK3)');

[err,act_f,act_frs] = asgard(fokkerplanck1_4p1b,'lev',4,'quiet',true,'deg',3);

verifyLessThan(testCase,err,7e-4);

end

function asgard_fokkerplanck4p1b_implicit_test(testCase)

addpath(genpath(pwd));

disp('Testing fokkerplanck4p1b (CN)');

[err,act_f,act_frs] = asgard(fokkerplanck1_4p1b,...
    'lev',4,'quiet',true,'deg',3,'timestep_method', 'CN');

verifyLessThan(testCase,err,7e-4);

end

function asgard_fokkerplanck4p1b_adapt_test(testCase)

addpath(genpath(pwd));

disp('Testing fokkerplanck4p1b (CN/adapt)');

[err,act_f,act_frs] = asgard(fokkerplanck1_4p1b,...
    'lev',4,'quiet',true,'deg',3,'timestep_method', 'CN','adapt',true);

verifyLessThan(testCase,err,1.5e-3); % TODO : need to look into why this is larger than the adapt=false approach.

end

function asgard_diffusion2_explicit_test(testCase)

addpath(genpath(pwd));

disp('Testing diffusion2 (RK3)');

[err,act_f,act_frs] = asgard(diffusion2,'lev',3,'quiet',true,'deg',3);

verifyLessThan(testCase,err,1.5e-5);

end

function asgard_diffusion2_implicit_test(testCase)

addpath(genpath(pwd));

disp('Testing diffusion2 (CN)');

[err,act_f,act_frs] = asgard(diffusion2, ...
    'lev',3,'quiet',true,'deg',3,'timestep_method', 'CN');

verifyLessThan(testCase,err,1.5e-5);

end

function asgard_diffusion2_oldhash_test(testCase)

addpath(genpath(pwd));

disp('Testing diffusion2 (CN/oldhash)');

[err,act_f,act_frs] = asgard(diffusion2, ...
    'lev',3,'quiet',true,'deg',3,'timestep_method', 'CN','use_oldhash',true);

verifyLessThan(testCase,err,1.5e-5);

end

function asgard_diffusion2_adapt_test(testCase)

addpath(genpath(pwd));

disp('Testing diffusion2 (CN/adapt)');

[err,act_f,act_frs] = asgard(diffusion2,'lev',3,'quiet',true,'deg',3, ...
    'timestep_method', 'CN','adapt',true);

verifyLessThan(testCase,err,9.5e-4);

end

function asgard_fokkerplanck1_5p1a_noLHS_implicit_test(testCase)

addpath(genpath(pwd));

disp('Testing fokkerplanck1_5p1a_noLHS (CN)');

[err,act_f,act_frs] = asgard(fokkerplanck1_5p1a_noLHS,...
    'timestep_method', 'CN','lev',3,'num_steps',30,'CFL',1.5,'quiet',true);

verifyLessThan(testCase,err,2.1e-2);

end

function asgard_fokkerplanck1_5p1a_implicit_test(testCase)

addpath(genpath(pwd));

disp('Testing fokkerplanck1_5p1a (BE / with LHS)');

[err,act_f,act_frs] = asgard(fokkerplanck1_5p1a,'timestep_method','BE', ...
    'lev',3,'num_steps',30,'CFL',1.5,'quiet',true);

verifyLessThan(testCase,err,2.5e-2);

end

