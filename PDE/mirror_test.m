function tests = mirror_test()

fh = localfunctions;
tests = functiontests(fh);

end

function mirror1_velocity_test(testCase)
addpath(genpath(pwd));
disp('Testing the velocity dimension within mirror1');
% setup PDE
args = {'lev',4,'deg',4,'dt',1e-8,'quiet',true,'num_steps',20,'timestep_method','matrix_exponential','normalize_by_mass',true};
opts = OPTS(args);
pde = mirror1_velocity(opts);
% modify PDE - not needed here
% run PDE
[err,fval,fval_realspace,nodes,err_realspace] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = err / norm(fval);
verifyLessThan(testCase,rel_err,1e-5);
end

function mirror1_pitch_test(testCase)
addpath(genpath(pwd));
disp('Testing the pitch dimension within mirror1');
% setup PDE
args = {'lev',4,'deg',4,'dt',5e-7,'quiet',true,'num_steps',20,'timestep_method','matrix_exponential','normalize_by_mass',false};
opts = OPTS(args);
pde = mirror1_pitch(opts);
% modify PDE - not needed here
% run PDE
[err,fval,fval_realspace,nodes,err_realspace] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = err / norm(fval);
verifyLessThan(testCase,rel_err,1e-5);
end

function mirror3_pitch_test(testCase)
addpath(genpath(pwd));
disp('Testing the pitch dimension within mirror3');
% setup PDE
args = {'lev',3,'deg',3,'dt',1e-10,'calculate_mass', true,'quiet',false,'num_steps',5,'normalize_by_mass', true, 'timestep_method','matrix_exponential'};
opts = OPTS(args);
pde = mirror3(opts);
% modify PDE
pde.dimensions{1}.init_cond_fn = @(x,p,t) x.*0 + 1;
pde.dimensions{2}.init_cond_fn = @(z,p,t) cos(z);
pde.dimensions{3}.init_cond_fn = @(x,p,t) x.*0 + 1;
pde.params.boundary_cond_v = @(v,p,t) v.*0 + 1;
pde.params.boundary_cond_s = @(s,p,t) s.*0 + 1;
pde.analytic_solutions_1D = { ...    
    @(v,p,t) exp(-p.nu_D(v,p.a,p.b).*t), ...
    @(z,p,t) cos(z), ...
    @(s,p,t) s.*0 + 1, ...
    @(t,p) t.*0 + 1; %pitch_t(t)
    };
% run PDE
[err,fval,fval_realspace,nodes,err_realspace] = asgard_run_pde(opts,pde);
% assert on correctness
verifyLessThan(testCase,err,1e-5);
end

function mirror3_velocity_test(testCase)

addpath(genpath(pwd));
disp('Testing the velocity dimension within mirror3');
% setup PDE
args = {'lev',3,'deg',5,'dt',1e-8, 'calculate_mass', true, 'normalize_by_mass', true, 'timestep_method', 'matrix_exponential','quiet',true,'num_steps',20};
opts = OPTS(args);
pde = mirror3(opts);
% modify PDE
pde.dimensions{1}.init_cond_fn = @(x,p,t) p.init_cond_v(x);
pde.dimensions{2}.init_cond_fn = @(x,p,t) x.*0 + 1;
pde.dimensions{3}.init_cond_fn = @(x,p,t) x.*0 + 1;
pde.params.boundary_cond_v = @(v,p,t) p.init_cond_v(v);
pde.params.boundary_cond_z = @(z,p,t) z.*0 + 1;
pde.params.boundary_cond_s = @(s,p,t) s.*0 + 1;
pde.dimensions{1,1}.min = 0;
pde.dimensions{1,1}.max = 5e6;
pde.dimensions{1,2}.min = 0;
pde.dimensions{1,2}.max = pi/2;

    function ret = solution(v,p,t)
        ret =  p.analytic_solution_v(v,p,t);
        if isfield(p,'norm_fac')
            ret = p.norm_fac .* ret;
        end
    end

pde.analytic_solutions_1D = { ...    
    @(v,p,t) solution(v,p,t), ...
    @(z,p,t) z.*0 + 1, ...
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
args = {'lev',3,'deg',4,'dt',1e-15, 'normalize_by_mass', true, 'timestep_method', 'matrix_exponential','quiet',true,'num_steps',5};
opts = OPTS(args);
pde = mirror3(opts);
% modify PDE
pde.dimensions{1}.init_cond_fn = @(x,p,t) x.*0 + 1;
pde.dimensions{2}.init_cond_fn = @(x,p,t) x.*0 + 1;
pde.dimensions{3}.init_cond_fn = @(x,p,t) exp(x);

pde.params.boundary_cond_v = @(v,p,t) v.*0 + 1;
pde.params.boundary_cond_z = @(z,p,t) z.*0 + 1;

    function ret = solution(s,p,t)
        ret =  p.analytic_solution_s(s,p,t);
        if isfield(p,'norm_fac')
            ret = p.norm_fac .* ret;
        end
    end

pde.analytic_solutions_1D = { ...    
    @(v,p,t) v.*0 + 1, ...
    @(z,p,t) z.*0 + 1, ...
    @(s,p,t) solution(s,p,t), ...
    @(t,p) t.*0 + 1; %pitch_t(t)
    };
% run PDE
[err,fval,fval_realspace,nodes,err_realspace] = asgard_run_pde(opts,pde);
% assert on correctness
verifyLessThan(testCase,err,1e-2);
end
