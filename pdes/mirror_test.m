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
args = {'lev',3,'deg',3,'dt',1e-7,'calculate_mass',true,'quiet',true,'num_steps',1,'normalize_by_mass',true,'timestep_method','matrix_exponential'};
opts = OPTS(args);
pde = mirror3(opts);
num_dims = numel(pde.dimensions);

% modify PDE
pde.params.boundary_cond_v = @(v,p,t) v.*0 + 1;
pde.params.boundary_cond_s = @(s,p,t) s.*0 + 1;
pde.terms = {pde.terms{1,3}};

s_v = @(v,p,t) exp(-p.nu_D(v,p.a,p.b).*t);
s_z = @(z,p,t) cos(z);
soln1 = new_md_func(num_dims,{s_v,s_z,[]});
pde.solutions = {soln1};

ic_z = @(z,p,t) cos(z);
ic1 = new_md_func(num_dims,{[],ic_z,[]});
pde.initial_conditions = {ic1};

% run PDE
[err,fval,fval_realspace,nodes,err_realspace,outputs] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = outputs.rel_err{end};
verifyLessThan(testCase,rel_err,1e-4);
end

function mirror3_velocity_test(testCase)

addpath(genpath(pwd));
disp('Testing the velocity dimension within mirror3');
% setup PDE
args = {'lev',3,'deg',3,'dt',10e-8,'calculate_mass',true,'normalize_by_mass',true,'timestep_method','matrix_exponential','quiet',true,'num_steps',1};
opts = OPTS(args);
pde = mirror3(opts);
num_dims = numel(pde.dimensions);

% modify PDE
ic_v = @(x,p,t) p.init_cond_v(x);
ic1 = new_md_func(num_dims,{ic_v,[],[]});
pde.initial_conditions = {ic1};

pde.params.boundary_cond_v = @(v,p,t) p.init_cond_v(v);
pde.params.boundary_cond_z = @(z,p,t) z.*0 + 1;
pde.params.boundary_cond_s = @(s,p,t) s.*0 + 1;

pde.dimensions{1,1}.min = 0;
pde.dimensions{1,1}.max = 5e6;
pde.dimensions{1,2}.min = 0;
pde.dimensions{1,2}.max = pi/2;
pde.terms = {pde.terms{1,1}, pde.terms{1,2}};

    function ret = solution(v,p,t)
        ret =  p.soln_v(v,p,t);
        if isfield(p,'norm_fac')
            ret = p.norm_fac .* ret;
        end
    end

soln1 = new_md_func(num_dims,{@solution,[],[]});
pde.solutions = {soln1};

% run PDE
[err,fval,fval_realspace,nodes,err_realspace,outputs] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = outputs.rel_err{end};
verifyLessThan(testCase,rel_err,1.5e-4);
end

function mirror3_space_test(testCase)
addpath(genpath(pwd));
disp('Testing the spatial dimension within mirror3');
% setup PDE
args = {'lev',3,'deg',3,'dt',5e-15, 'normalize_by_mass', true, 'timestep_method', 'matrix_exponential','quiet',true,'num_steps',1};
opts = OPTS(args);
pde = mirror3(opts);
num_dims = numel(pde.dimensions);

% modify PDE
ic_s = @(x,p,t) exp(x);
ic1 = new_md_func(num_dims,{[],[],ic_s});
pde.initial_conditions = {ic1};

pde.params.boundary_cond_v = @(v,p,t) v.*0 + 1;
pde.params.boundary_cond_z = @(z,p,t) z.*0 + 1;
pde.terms = {pde.terms{1,4}, pde.terms{1,5}};

    function ret = solution(s,p,t)
        ret =  p.soln_s(s,p,t);
        if isfield(p,'norm_fac')
            ret = p.norm_fac .* ret;
        end
    end

soln1 = new_md_func(num_dims,{[],[],@solution});
pde.solutions = {soln1};

% run PDE
[err,fval,fval_realspace,nodes,err_realspace,outputs] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = outputs.rel_err{end};
verifyLessThan(testCase,rel_err,1e-4);
end

function mirror3_manufacture_test(testCase)
addpath(genpath(pwd));
disp('Testing a manufactured solution within mirror3');
% setup PDE
args = {'lev',3,'deg',3,'dt',1e-12,'calculate_mass',true,'quiet',true,'num_steps',1,'normalize_by_mass',true,'timestep_method','BE'};
opts = OPTS(args);
pde = mirror3(opts);
num_dims = numel(pde.dimensions);

% modify PDE
pde.params.nu_D = @(v,a,b) v;
pde.params.nu_s = @(v,a,b) v;
pde.params.nu_par = @(v,a,b) v;
pde.params.boundary_cond_v = @(v,p,t) exp(-v.*t);
pde.params.boundary_cond_s = @(s,p,t) s.*0 + 1;
pde.terms = {pde.terms{1,1},pde.terms{1,2},pde.terms{1,3},pde.terms{1,4}};

%source function for manufactured solution
src_v = @(v,p,t) exp(-v.*t).*((p.a.m./(p.a.m + p.b.m)).*(-4*v + v.^2.*t) + 5*v.^2.*t./2 - v.^3.*t.^2./2);
src_z = @(z,p,t) cos(z);

s_v = @(v,p,t) exp(-v.*t);
s_z = @(z,p,t) cos(z);
soln1 = new_md_func(num_dims,{s_v,s_z,[]});
src = new_md_func(num_dims,{src_v,src_z,[]});
pde.solutions = {soln1};
pde.sources = {src};

ic_z = @(z,p,t) cos(z);
ic1 = new_md_func(num_dims,{[],ic_z,[]});
pde.initial_conditions = {ic1};

% run PDE
[err,fval,fval_realspace,nodes,err_realspace,outputs] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = outputs.rel_err{end};
verifyLessThan(testCase,rel_err,1e-4);
end
