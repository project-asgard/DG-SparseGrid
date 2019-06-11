function tests = asgard_test()

fh = localfunctions;
tests = functiontests(fh);

end

function asgard_continuity1_explicit_test(testCase)

addpath(genpath(pwd));

disp('Testing continuity1 (explicit)');

[err,act_f,act_frs] = asgard(continuity1,'quiet',true,'num_steps',2);

verifyLessThan(testCase,err,1e-4);

end

function asgard_continuity1_implicit_test(testCase)

addpath(genpath(pwd));

disp('Testing continuity1 (implicit)');

[err,act_f,act_frs] = asgard(continuity1,'quiet',true,'implicit',true,'num_steps',2);

verifyLessThan(testCase,err,1e-4);

end

function asgard_continuity1_adapt_test(testCase)

addpath(genpath(pwd));

disp('Testing continuity1 (adapt)');

[err,act_f,act_frs] = asgard(continuity1,'quiet',true,'implicit',true,'num_steps',2,'adapt',true);

verifyLessThan(testCase,err,1e-4);

end

function asgard_continuity1_oldhash_test(testCase)

addpath(genpath(pwd));

disp('Testing continuity1 (old_hash)');

[err,act_f,act_frs] = asgard(continuity1,'quiet',true,'num_steps',2,'use_oldhash',true);

verifyLessThan(testCase,err,1e-4);

end

function asgard_continuity2_expliciti_test(testCase)

addpath(genpath(pwd));

disp('Testing continuity2 (explicit)');

[err,act_f,act_frs] = asgard(continuity2,'quiet',true,'lev',4,'deg',3,'num_steps',2);

verifyLessThan(testCase,err,1e-4);

end

function asgard_continuity2_implicit_test(testCase)

addpath(genpath(pwd));

disp('Testing continuity2 (implicit)');

[err,act_f,act_frs] = asgard(continuity2,'quiet',true,'implicit',true,'lev',4,'deg',3,'num_steps',2);

verifyLessThan(testCase,err,1e-4);

end

function asgard_continuity2_adapt_test(testCase)

addpath(genpath(pwd));

disp('Testing continuity2 (adapt)');

[err,act_f,act_frs] = asgard(continuity2,'quiet',true,'implicit',true,'lev',4,'deg',3,'num_steps',2,'adapt',true);

verifyLessThan(testCase,err,1e-4);

end


