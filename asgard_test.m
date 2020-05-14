function tests = asgard_test()

fh = localfunctions;
tests = functiontests(fh);

end

function asgard_advection1_explicit_test(testCase)

addpath(genpath(pwd));

disp('Testing advection1 (explicit)');

[err, act_f, act_frs] = asgard(advection1, 'quiet', true, 'num_steps', 5);

verifyLessThan(testCase, err, 9e-4);

end

function asgard_advection1reverse_explicit_test(testCase)

addpath(genpath(pwd));

disp('Testing advection1_reverse_flow (explicit)');

[err, act_f, act_frs] = asgard(advection1_reverse_flow, 'quiet', true, 'num_steps', 5);

verifyLessThan(testCase, err, 9e-4);

end

function asgard_advection1_implicit_test(testCase)

addpath(genpath(pwd));

disp('Testing advection1 (implicit)');

[err, act_f, act_frs] = asgard(advection1, 'quiet',true,'implicit',true,'lev',4,'deg',3,'num_steps',5);

verifyLessThan(testCase, err, 3e-6);

end

function asgard_advection1reverse_implicit_test(testCase)

addpath(genpath(pwd));

disp('Testing advection1_reverse_flow (implicit)');

[err, act_f, act_frs] = asgard(advection1_reverse_flow, 'quiet',true,'implicit',true,'lev',4,'deg',3,'num_steps',5);

verifyLessThan(testCase, err, 3e-6);

end

function asgard_diffusion1_explicit_test(testCase)

addpath(genpath(pwd));

disp('Testing diffusion1 (explicit)');

[err,act_f,act_frs] = asgard(diffusion1,'lev',3,'quiet',true,'deg',3);

verifyLessThan(testCase,err,2.8e-5);

end

function asgard_diffusion1_implicit_test(testCase)

addpath(genpath(pwd));

disp('Testing diffusion1 (implicit)');

[err,act_f,act_frs] = asgard(diffusion1,'lev',3,'quiet',true,'deg',3,'implicit',true);

verifyLessThan(testCase,err,2.8e-5);

end


function asgard_continuity1_ex_test(testCase)

addpath(genpath(pwd));

disp('Testing continuity1 (explicit)');

[err,act_f,act_frs] = asgard(continuity1,'quiet',true,'num_steps',2);

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

function asgard_fokkerplanck4p1b_explicit_test(testCase)

addpath(genpath(pwd));

disp('Testing fokkerplanck4p1b (explicit)');

[err,act_f,act_frs] = asgard(fokkerplanck1_4p1b,'lev',4,'quiet',true,'deg',3);

verifyLessThan(testCase,err,7e-4);

end

function asgard_fokkerplanck4p1b_implicit_test(testCase)

addpath(genpath(pwd));

disp('Testing fokkerplanck4p1b (implicit)');

[err,act_f,act_frs] = asgard(fokkerplanck1_4p1b,'lev',4,'quiet',true,'deg',3,'implicit',true);

verifyLessThan(testCase,err,7e-4);

end

function asgard_fokkerplanck4p1b_adapt_test(testCase)

addpath(genpath(pwd));

disp('Testing fokkerplanck4p1b (adapt)');

[err,act_f,act_frs] = asgard(fokkerplanck1_4p1b,'lev',4,'quiet',true,'deg',3,'implicit',true,'adapt',true);

verifyLessThan(testCase,err,1.5e-3); % TODO : need to look into why this is larger than the adapt=false approach.

end

function asgard_diffusion2_explicit_test(testCase)

addpath(genpath(pwd));

disp('Testing diffusion2 (explicit)');

[err,act_f,act_frs] = asgard(diffusion2,'lev',3,'quiet',true,'deg',3);

verifyLessThan(testCase,err,1.5e-5);

end

function asgard_diffusion2_implicit_test(testCase)

addpath(genpath(pwd));

disp('Testing diffusion2 (implicit)');

[err,act_f,act_frs] = asgard(diffusion2,'lev',3,'quiet',true,'deg',3,'implicit',true);

verifyLessThan(testCase,err,1.5e-5);

end

function asgard_diffusion2_oldhash_test(testCase)

addpath(genpath(pwd));

disp('Testing diffusion2 (oldhash)');

[err,act_f,act_frs] = asgard(diffusion2,'lev',3,'quiet',true,'deg',3,'implicit',true,'use_oldhash',true);

verifyLessThan(testCase,err,1.5e-5);

end

function asgard_diffusion2_adapt_test(testCase)

addpath(genpath(pwd));

disp('Testing diffusion2 (adapt)');

[err,act_f,act_frs] = asgard(diffusion2,'lev',3,'quiet',true,'deg',3,'implicit',true,'adapt',true);

verifyLessThan(testCase,err,9.5e-4);

end

function asgard_fokkerplanck1_5p1a_noLHS_implicit_test(testCase)

addpath(genpath(pwd));

disp('Testing fokkerplanck1_5p1a_noLHS (implicit)');

[err,act_f,act_frs] = asgard(fokkerplanck1_5p1a_noLHS,'implicit',true,'lev',3,'num_steps',30,'CFL',1.5,'quiet',true);

verifyLessThan(testCase,err,2.1e-2);

end

function asgard_fokkerplanck1_5p1a_implicit_test(testCase)

addpath(genpath(pwd));

disp('Testing fokkerplanck1_5p1a (implicit / with LHS)');

[err,act_f,act_frs] = asgard(fokkerplanck1_5p1a,'implicit',true,'lev',3,'num_steps',30,'CFL',1.5,'quiet',true,'implicit_method','BE');

verifyLessThan(testCase,err,2.5e-2);

end


function asgard_advection1_implicit_BDF2_test(testCase)

addpath(genpath(pwd));

disp('Testing advection1 (implicit) BDF2');

[err, act_f, act_frs] = asgard(fokkerplanck1_4p1a, 'quiet',true,'implicit',true,'lev',5,'deg',3,'num_steps',1, 'dt', 0.01, 'implicit_method', 'BDF2');

[err1, acf_f1, act_frs1] = asgard(fokkerplanck1_4p1a, 'quiet',true,'implicit',true,'lev',5,'deg',3,'num_steps',2, 'dt', 0.005, 'implicit_method', 'BDF2');

[err2, act_f2, act_frs2] = asgard(fokkerplanck1_4p1a, 'quiet',true,'implicit',true,'lev',5,'deg',3,'num_steps',4, 'dt', 0.0025, 'implicit_method', 'BDF2');

logslope1 = (log10(err1) - log10(err))/(log10(1/0.01) - log10(1/0.005));

logslope2 = (log10(err2) - log10(err1))/(log10(1/0.005) - log10(1/0.0025));

%slope_err1 = abs(logslope1 - 2.000);
slope_err2 = abs(logslope2 - 2.000);

%verifyLessThan(testCase,slope_err1, 1e-5);

verifyLessThan(testCase,slope_err2, 1e-6);


end


