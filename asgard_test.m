function tests = asgard_test()

fh = localfunctions;
tests = functiontests(fh);

end

function asgard_advection1_explicit_test(testCase)
addpath(genpath(pwd));
disp('Testing advection1 (RK3)');
[err, act_f, act_frs] = asgard(@advection1,'quiet',true,'num_steps',5);
verifyLessThan(testCase, err, 9e-4);
end

function asgard_advection1reverse_explicit_test(testCase)
addpath(genpath(pwd));
disp('Testing advection1_reverse_flow (RK3)');
[err, act_f, act_frs] = asgard(@advection1_reverse_flow,'quiet',true,'num_steps', 5);
verifyLessThan(testCase, err, 9e-4);
end

function asgard_advection1_implicit_test(testCase)
addpath(genpath(pwd));
disp('Testing advection1 (CN)');
[err, act_f, act_frs] = asgard(@advection1,'quiet',true,...
    'timestep_method','CN','lev',4,'deg',3,'num_steps',5);
verifyLessThan(testCase, err, 3e-6);
end

function asgard_advection1_implicit_BE_test(testCase)
addpath(genpath(pwd));
disp('Testing advection1 (BE)');
[err, act_f, act_frs] = asgard(@advection1,'quiet',true,...
    'timestep_method','BE','lev',4,'deg',3,'num_steps',5);
verifyLessThan(testCase, err, 3e-6);
end

function asgard_advection1_ode45_test(testCase)
addpath(genpath(pwd));
disp('Testing advection1 (ode45)');
[err, act_f, act_frs] = asgard(@advection1,'quiet',true,...
    'timestep_method','ode45','lev',4,'deg',3,'num_steps',1,'dt',0.0019635*5);
verifyLessThan(testCase, err, 3e-6);
end

function asgard_advection1_ode15s_test(testCase)
addpath(genpath(pwd));
disp('Testing advection1 (ode15s)');
[err, act_f, act_frs] = asgard(@advection1,'quiet',true,...
    'timestep_method','ode15s','lev',4,'deg',3,'num_steps',1,'dt',0.0019635*5);
verifyLessThan(testCase, err, 3e-6);
end

function asgard_advection1reverse_implicit_test(testCase)
addpath(genpath(pwd));
disp('Testing advection1_reverse_flow (CN)');
[err, act_f, act_frs] = asgard(@advection1_reverse_flow, ...
    'quiet',true,'timestep_method', 'CN','lev',4,'deg',3,'num_steps',5);
verifyLessThan(testCase, err, 3e-6);
end

function asgard_advection1reverse_implicit_BE_test(testCase)
addpath(genpath(pwd));
disp('Testing advection1_reverse_flow (BE)');
[err, act_f, act_frs] = asgard(@advection1_reverse_flow, ...
    'quiet',true,'timestep_method', 'BE','lev',4,'deg',3,'num_steps',5);
verifyLessThan(testCase, err, 3e-6);
end

function asgard_diffusion1_explicit_test(testCase)
addpath(genpath(pwd));
disp('Testing diffusion1 (RK3)');
[err,act_f,act_frs] = asgard(@diffusion1,'lev',3,'quiet',true,'deg',3);
verifyLessThan(testCase,err,2.8e-5);
end

function asgard_diffusion1_ode45_test(testCase)
addpath(genpath(pwd));
disp('Testing diffusion1 (ode45)');
[err,act_f,act_frs] = asgard(@diffusion1,'lev',3,'quiet',true,...
    'deg',3,'dt',0.00015625*5,'num_steps',1,'timestep_method','ode45');
verifyLessThan(testCase,err,2.8e-5);
end

function asgard_diffusion1_ode15s_test(testCase)
addpath(genpath(pwd));
disp('Testing diffusion1 (ode15s)');
[err,act_f,act_frs] = asgard(@diffusion1,'lev',3,'quiet',true,...
    'deg',3,'dt',0.00015625*5,'num_steps',1,'timestep_method','ode15s');
verifyLessThan(testCase,err,2.8e-5);
end

function asgard_diffusion1_ode15i_test(testCase)
addpath(genpath(pwd));
disp('Testing diffusion1 (ode15i)');
[err,act_f,act_frs] = asgard(@diffusion1,'lev',3,'quiet',true,...
    'deg',3,'dt',0.00015625*5,'num_steps',1,'timestep_method','ode15i');
verifyLessThan(testCase,err,2.8e-5);
end

function asgard_diffusion1_implicit_test(testCase)
addpath(genpath(pwd));
disp('Testing diffusion1 (CN)');
[err,act_f,act_frs] = asgard(@diffusion1,'lev',3,'quiet',true,'deg',3,'timestep_method', 'CN');
verifyLessThan(testCase,err,2.8e-5);
end

function asgard_diffusion1_implicit_BE_test(testCase)
addpath(genpath(pwd));
disp('Testing diffusion1 (BE)');
[err,act_f,act_frs] = asgard(@diffusion1,'lev',3,'quiet',true,'deg',3,'timestep_method', 'BE');
verifyLessThan(testCase,err,2.8e-5);
end

function asgard_continuity1_ex_test(testCase)
addpath(genpath(pwd));
disp('Testing continuity1 (RK3)');
[err,act_f,act_frs] = asgard(@continuity1,'quiet',true,'num_steps',2);
verifyLessThan(testCase,err,1e-4);
end


function asgard_continuity1_adapt_test(testCase)
addpath(genpath(pwd));
disp('Testing continuity1 (CN/adapt)');
[err,act_f,act_frs] = asgard(@continuity1,'quiet',true,...
    'timestep_method', 'CN','num_steps',2,'adapt',true);
verifyLessThan(testCase,err,1e-4);
end

function asgard_continuity1_oldhash_test(testCase)
addpath(genpath(pwd));
disp('Testing continuity1 (RK3/old_hash)');
[err,act_f,act_frs] = asgard(@continuity1,'quiet',true,'num_steps',2,'use_oldhash',true);
verifyLessThan(testCase,err,1e-4);
end

function asgard_continuity2_expliciti_test(testCase)
addpath(genpath(pwd));
disp('Testing continuity2 (RK3)');
[err,act_f,act_frs] = asgard(@continuity2,'quiet',true,'lev',4,'deg',3,'num_steps',2);
verifyLessThan(testCase,err,1e-4);
end

function asgard_continuity2_implicit_CN_test(testCase)
addpath(genpath(pwd));
disp('Testing continuity2 (CN)');
[err,act_f,act_frs] = asgard(@continuity2,'quiet',true,...
    'timestep_method', 'CN','lev',4,'deg',3,'num_steps',2);
verifyLessThan(testCase,err,1e-4);
end

function asgard_continuity2_implicit_BE_test(testCase)
addpath(genpath(pwd));
disp('Testing continuity2 (BE)');
[err,act_f,act_frs] = asgard(@continuity2,'quiet',true,...
    'timestep_method', 'BE','lev',4,'deg',3,'num_steps',2);
verifyLessThan(testCase,err,1e-4);
end

function asgard_continuity2_adapt_test(testCase)
addpath(genpath(pwd));
disp('Testing continuity2 (CN/adapt)');
[err,act_f,act_frs] = asgard(@continuity2,'quiet',true,...
    'timestep_method', 'CN','lev',4,'deg',3,'num_steps',2,'adapt',true);
verifyLessThan(testCase,err,1e-4);
end

function asgard_fokkerplanck1_pitch_E_explicit_test(testCase)
addpath(genpath(pwd));
disp('Testing fokkerplanck1_pitch_E (RK3)');
[err,act_f,act_frs] = asgard(@fokkerplanck1_pitch_E,'lev',4,'quiet',true,'deg',3,'max_lev',8,'case',2);
verifyLessThan(testCase,err,7e-4);
end

function asgard_fokkerplanck1_pitch_E_implicit_CN_test(testCase)
addpath(genpath(pwd));
disp('Testing fokkerplanck1_pitch_E (CN)');
[err,act_f,act_frs] = asgard(@fokkerplanck1_pitch_E,...
    'lev',4,'quiet',true,'deg',3,'timestep_method','CN','max_lev',8,'case',2);
verifyLessThan(testCase,err,7e-4);
end

function asgard_fokkerplanck1_pitch_E_implicit_BE_test(testCase)
addpath(genpath(pwd));
disp('Testing fokkerplanck1_pitch_E (BE)');
[err,act_f,act_frs] = asgard(@fokkerplanck1_pitch_E,...
    'lev',4,'quiet',true,'deg',3,'timestep_method','BE','max_lev',8,'case',2);
verifyLessThan(testCase,err,7e-4);
end

function asgard_fokkerplanck1_pitch_E_adapt_test(testCase)
addpath(genpath(pwd));
disp('Testing fokkerplanck1_pitch_E (CN/adapt)');
[err,act_f,act_frs] = asgard(@fokkerplanck1_pitch_E,...
    'lev',4,'quiet',true,'deg',3,'timestep_method','CN','adapt',true,'max_lev',8,'case',2);
verifyLessThan(testCase,err,1.5e-3); % TODO : need to look into why this is larger than the adapt=false approach.
end

function asgard_diffusion2_explicit_test(testCase)
addpath(genpath(pwd));
disp('Testing diffusion2 (RK3)');
[err,act_f,act_frs] = asgard(@diffusion2,'lev',3,'quiet',true,'deg',3);
verifyLessThan(testCase,err,1.5e-5);
end

function asgard_diffusion2_implicit_CN_test(testCase)
addpath(genpath(pwd));
disp('Testing diffusion2 (CN)');
[err,act_f,act_frs] = asgard(@diffusion2, ...
    'lev',3,'quiet',true,'deg',3,'timestep_method', 'CN');
verifyLessThan(testCase,err,1.5e-5);
end

function asgard_diffusion2_implicit_BE_test(testCase)
addpath(genpath(pwd));
disp('Testing diffusion2 (BE)');
[err,act_f,act_frs] = asgard(@diffusion2, ...
    'lev',3,'quiet',true,'deg',3,'timestep_method', 'BE');
verifyLessThan(testCase,err,1.5e-5);
end

function asgard_diffusion2_oldhash_test(testCase)
addpath(genpath(pwd));
disp('Testing diffusion2 (CN/oldhash)');
[err,act_f,act_frs] = asgard(@diffusion2, ...
    'lev',3,'quiet',true,'deg',3,'timestep_method', 'CN','use_oldhash',true);
verifyLessThan(testCase,err,1.5e-5);
end

function asgard_diffusion2_oldhash_connectivity_SG_test(testCase)
addpath(genpath(pwd));
disp('Testing diffusion2 (CN/oldhash/connectivity/SG)');
[err,act_f,act_frs] = asgard(@diffusion2,'lev',3,'quiet',true,'deg',3,'timestep_method','BE','use_oldhash',true,'use_connectivity',true);
verifyLessThan(testCase,err,1.4e-5);
end

function asgard_diffusion2_oldhash_connectivity_FG_test(testCase)
addpath(genpath(pwd));
disp('Testing diffusion2 (CN/oldhash/connectivity/FG)');
[err,act_f,act_frs] = asgard(@diffusion2,'lev',3,'quiet',true,'deg',3,'timestep_method','BE','use_oldhash',true,'use_connectivity',true,'grid_type','FG');
verifyLessThan(testCase,err,1.4e-5);
end

function asgard_diffusion2_adapt_test(testCase)
addpath(genpath(pwd));
disp('Testing diffusion2 (CN/adapt)');
[err,act_f,act_frs] = asgard(@diffusion2,'lev',3,'quiet',true,'deg',3, ...
    'timestep_method', 'CN','adapt',true);
verifyLessThan(testCase,err,9.5e-4);
end

function asgard_fokkerplanck1_5p1a_noLHS_ode15s_test(testCase)
addpath(genpath(pwd));
disp('Testing fokkerplanck1_5p1a_noLHS (ode15s)');
[err,act_f,act_frs] = asgard(@fokkerplanck1_5p1a_noLHS,...
    'timestep_method','ode15s','lev',3,'num_steps',1,'dt',1.875*30,'quiet',true,'case',2);
verifyLessThan(testCase,err,3.8e-2);
end

function asgard_fokkerplanck1_5p1a_noLHS_CN_test(testCase)
addpath(genpath(pwd));
disp('Testing fokkerplanck1_5p1a_noLHS (CN)');
[err,act_f,act_frs] = asgard(@fokkerplanck1_5p1a_noLHS,...
    'timestep_method','CN','lev',3,'num_steps',30,'dt',1.875,'quiet',true,'case',2);
verifyLessThan(testCase,err,3.8e-2);
end

function asgard_fokkerplanck1_5p1a_noLHS_CN_TIA_test(testCase)
addpath(genpath(pwd));
disp('Testing fokkerplanck1_5p1a_noLHS (CN / TIA)');
[err,act_f,act_frs] = asgard(@fokkerplanck1_5p1a_noLHS,...
    'timestep_method','CN','lev',3,'num_steps',30,'dt',1.875,...
    'quiet',true,'time_independent_A',true,'case',2);
verifyLessThan(testCase,err,3.8e-2);
end

function asgard_fokkerplanck1_5p1a_noLHS_BE_test(testCase)
addpath(genpath(pwd));
disp('Testing fokkerplanck1_5p1a_noLHS (BE)');
[err,act_f,act_frs] = asgard(@fokkerplanck1_5p1a_noLHS,...
    'timestep_method','BE','lev',3,'num_steps',30,'dt',1.875,'quiet',true,'case',2);
verifyLessThan(testCase,err,3.8e-2);
end

function asgard_fokkerplanck1_5p1a_noLHS_BE_TIA_test(testCase)
addpath(genpath(pwd));
disp('Testing fokkerplanck1_5p1a_noLHS (BE / TIA)');
[err,act_f,act_frs] = asgard(@fokkerplanck1_5p1a_noLHS,...
    'timestep_method','BE','lev',3,'num_steps',30,'dt',1.875,...
    'quiet',true,'time_independent_A',true,'case',2);
verifyLessThan(testCase,err,3.8e-2);
end

function asgard_fokkerplanck1_5p1a_implicit_BE_test(testCase)
addpath(genpath(pwd));
disp('Testing fokkerplanck1_5p1a (BE / with LHS)');
[err,act_f,act_frs] = asgard(@fokkerplanck1_5p1a,'timestep_method','BE', ...
    'lev',3,'num_steps',30,'CFL',1.5,'quiet',true);
verifyLessThan(testCase,err,2.5e-2);
end

function asgard_fokkerplanck1_5p1a_implicit_BE_TIA_test(testCase)
addpath(genpath(pwd));
disp('Testing fokkerplanck1_5p1a (BE / with LHS / TIA)');
[err,act_f,act_frs] = asgard(@fokkerplanck1_5p1a,'timestep_method','BE', ...
    'lev',3,'num_steps',30,'CFL',1.5,'quiet',true,'time_independent_A',true);
verifyLessThan(testCase,err,2.5e-2);
end

function asgard_fokkerplanck1_5p1a_implicit_CN_test(testCase)
addpath(genpath(pwd));
disp('Testing fokkerplanck1_5p1a (CN / with LHS)');
[err,act_f,act_frs] = asgard(@fokkerplanck1_5p1a,'timestep_method','CN', ...
    'lev',3,'num_steps',30,'CFL',1.5,'quiet',true);
verifyLessThan(testCase,err,2.5e-2);
end

function asgard_fokkerplanck1_5p1a_implicit_CN_TIA_test(testCase)
addpath(genpath(pwd));
disp('Testing fokkerplanck1_5p1a (CN / with LHS / TIA)');
[err,act_f,act_frs] = asgard(@fokkerplanck1_5p1a,'timestep_method','CN', ...
    'lev',3,'num_steps',30,'CFL',1.5,'quiet',true,'time_independent_A',true);
verifyLessThan(testCase,err,2.5e-2);
end

function end_points_grid_test(testCase)
addpath(genpath(pwd));
disp('Testing quadrature_with_end_points');
[err1,~,~,~,err_rsq]=asgard(@projecti_diff1,'timestep_method','BE', 'lev',3,'deg',3,'num_steps',1,'dt',0.0001,'quiet',true,'output_grid','quadrature');
[err2,~,~,~,err_rsf]=asgard(@projecti_diff1,'timestep_method','BE', 'lev',3,'deg',3,'num_steps',1,'dt',0.0001,'quiet',true,'output_grid','quadrature_with_end_points');
relErr = abs(err_rsq-err_rsf)/abs(max([err_rsq,err_rsf]));
verifyLessThan(testCase,relErr,1);
% dlg - i'm kind of surprised the difference between the quadrature points
% and those plus two end points give such a larger difference in the error,
% although it plots correctly.
end

function fixed_grid_test(testCase)
addpath(genpath(pwd));
disp('Testing fixed_grid routine');
[err1,~,~,~,err_rsq]=asgard(@projecti_diff1,'timestep_method','BE', 'lev',3,'deg',3,'num_steps',1,'dt',0.0001,'quiet',true,'output_grid','quadrature');
[err2,~,~,~,err_rsf]=asgard(@projecti_diff1,'timestep_method','BE', 'lev',3,'deg',3,'num_steps',1,'dt',0.0001,'quiet',true,'output_grid','fixed');
relErr = abs(err_rsq-err_rsf)/abs(max([err_rsq,err_rsf]));
verifyLessThan(testCase,relErr,1);
end

function uniform_grid_test(testCase)
addpath(genpath(pwd));
disp('Testing uniform_grid routine on projecti_diff1');
[err1,~,~,~,err_rsq]=asgard(@projecti_diff1,'timestep_method','BE', 'lev',3,'deg',3,'num_steps',1,'dt',0.0001,'quiet',true,'output_grid','quadrature');
[err2,~,~,~,err_rsu]=asgard(@projecti_diff1,'timestep_method','BE', 'lev',3,'deg',3,'num_steps',1,'dt',0.0001,'quiet',true,'output_grid','uniform');
relErr = abs(err_rsq-err_rsu)/abs(max([err_rsq,err_rsu]));
verifyLessThan(testCase,relErr,1);
end

function end_points_grid_2D_test(testCase)
addpath(genpath(pwd));
disp('Testing output on end points on diffusion2');
[err1,~,~,~,err_rsq]=asgard(@diffusion2,'timestep_method','BE', 'lev',3,'deg',3,'num_steps',1,'dt',0.0001,'quiet',true,'output_grid','quadrature');
[err2,~,~,~,err_rsf]=asgard(@diffusion2,'timestep_method','BE', 'lev',3,'deg',3,'num_steps',1,'dt',0.0001,'quiet',true,'output_grid','quadrature_with_end_points');
relErr = abs(err_rsq-err_rsf)/abs(max([err_rsq,err_rsf]));
verifyLessThan(testCase,relErr,1);
end

function fixed_grid_2D_test(testCase)
addpath(genpath(pwd));
disp('Testing fixed_grid routine on diffusion2');
[err1,~,~,~,err_rsq]=asgard(@diffusion2,'timestep_method','BE', 'lev',3,'deg',3,'num_steps',1,'dt',0.0001,'quiet',true,'output_grid','quadrature');
[err2,~,~,~,err_rsf]=asgard(@diffusion2,'timestep_method','BE', 'lev',3,'deg',3,'num_steps',1,'dt',0.0001,'quiet',true,'output_grid','fixed');
relErr = abs(err_rsq-err_rsf)/abs(max([err_rsq,err_rsf]));
verifyLessThan(testCase,relErr,1);
end

function uniform_grid_2D_test(testCase)
addpath(genpath(pwd));
disp('Testing fixed_grid routine');
[err1,~,~,~,err_rsq]=asgard(@diffusion2,'timestep_method','BE', 'lev',3,'deg',3,'num_steps',1,'dt',0.0001,'quiet',true,'output_grid','quadrature');
[err2,~,~,~,err_rsu]=asgard(@diffusion2,'timestep_method','BE', 'lev',3,'deg',3,'num_steps',1,'dt',0.0001,'quiet',true,'output_grid','uniform');
relErr = abs(err_rsq-err_rsu)/abs(max([err_rsq,err_rsu]));
verifyLessThan(testCase,relErr,1);
end

function mirror3_velocity_test(testCase)

addpath(genpath(pwd));
disp('Testing the velocity dimension within mirror3');
% setup PDE
args = {'lev',3,'deg',7,'dt',5e-8, 'normalize_by_mass', true, 'timestep_method', 'BE','quiet',true,'num_steps',10};
opts = OPTS(args);
pde = mirror3(opts);
% modify PDE
pde.dimensions{1}.init_cond_fn = @(x,p,t) p.maxwell(x,1e6,2e5);
pde.dimensions{2}.init_cond_fn = @(x,p,t) x.*0 + 1;
pde.dimensions{3}.init_cond_fn = @(x,p,t) x.*0 + 1;
pde.analytic_solutions_1D = { ...    
    @(v,p,t) p.analytic_solution_v(v,p,t), ...
    @(z,p,t) z.*0 + 1, ...
    @(s,p,t) s.*0 + 1, ...
    @(t,p) t.*0 + 1; %pitch_t(t)
    };
% run PDE
[err,fval,fval_realspace,nodes,err_realspace] = asgard_run_pde(opts,pde);
% assert on correctness
verifyLessThan(testCase,err,1e-2);

end

function mirror3_pitch_test(testCase)
addpath(genpath(pwd));
disp('Testing the pitch dimension within mirror3');
% setup PDE
args = {'lev',4,'deg',3,'dt',1e-14, 'normalize_by_mass', true, 'timestep_method', 'BE','quiet',true,'num_steps',5};
opts = OPTS(args);
pde = mirror3(opts);
% modify PDE
pde.dimensions{1}.init_cond_fn = @(x,p,t) x.*0 + 1;
pde.dimensions{2}.init_cond_fn = @(z,p,t) cos(z/2);
pde.dimensions{3}.init_cond_fn = @(x,p,t) x.*0 + 1;
pde.analytic_solutions_1D = { ...    
    @(v,p,t) v.*0 + 1, ...
    @(z,p,t) p.analytic_solution_z(z,p,t), ...
    @(s,p,t) s.*0 + 1, ...
    @(t,p) t.*0 + 1; %pitch_t(t)
    };
% run PDE
[err,fval,fval_realspace,nodes,err_realspace] = asgard_run_pde(opts,pde);
% assert on correctness
verifyLessThan(testCase,err,1e-2);
end

function mirror3_space_test(testCase)
addpath(genpath(pwd));
disp('Testing the spatial dimension within mirror3');
% setup PDE
args = {'lev',4,'deg',4,'dt',1e-15, 'normalize_by_mass', true, 'timestep_method', 'BE','quiet',true,'num_steps',5};
opts = OPTS(args);
pde = mirror3(opts);
% modify PDE
pde.dimensions{1}.init_cond_fn = @(x,p,t) x.*0 + 1;
pde.dimensions{2}.init_cond_fn = @(x,p,t) x.*0 + 1;
pde.dimensions{3}.init_cond_fn = @(x,p,t) exp(x);
pde.analytic_solutions_1D = { ...    
    @(v,p,t) v.*0 + 1, ...
    @(z,p,t) z.*0 + 1, ...
    @(s,p,t) p.analytic_solution_s(s,p,t), ...
    @(t,p) t.*0 + 1; %pitch_t(t)
    };
% run PDE
[err,fval,fval_realspace,nodes,err_realspace] = asgard_run_pde(opts,pde);
% assert on correctness
verifyLessThan(testCase,err,1e-2);
end

