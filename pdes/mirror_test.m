function tests = mirror_test()

fh = localfunctions;
tests = functiontests(fh);

end

function mirror1_pitch_div_test(testCase)
addpath(genpath(pwd));
disp('Testing the pitch angle dimension within mirror1');
% setup PDE
args = {'lev',5,'deg',4,'dt',5e-7,'quiet',true,'num_steps',10,'timestep_method','matrix_exponential','normalize_by_mass',true, 'calculate_mass', true};
opts = OPTS(args);
pde = mirror1_pitch_div(opts);
num_dims = numel(pde.dimensions);
% modify PDE
% run PDE
[err,fval,fval_realspace,nodes,err_realspace] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = err / norm(fval);
verifyLessThan(testCase,rel_err,1e-5);
end

function mirror1_collision_div_sameTemp_test(testCase)
addpath(genpath(pwd));
disp('Testing the collision operator within mirror1');
% setup PDE
args = {'lev',5,'deg',4,'dt',1e-3,'quiet',true,'num_steps',2,'timestep_method','matrix_exponential','normalize_by_mass',true, 'calculate_mass', true};
opts = OPTS(args);
pde = mirror1_collision_div(opts);
num_dims = numel(pde.dimensions);
% modify PDE
ic_v = @(v,p,t) p.soln_v(v,p,t);
ic1 = new_md_func(num_dims,{ic_v});
pde.initial_conditions = {ic1};
% run PDE
[err,fval,fval_realspace,nodes,err_realspace] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = err / norm(fval);
verifyLessThan(testCase,rel_err,1e-5);
end

function mirror1_collision_div_diffTemp_test(testCase)
addpath(genpath(pwd));
disp('Testing the collision operator within mirror1 by changing the temperature of the beam');
% setup PDE
args = {'lev',5,'deg',4,'dt',1e-3,'quiet',true,'num_steps',10,'timestep_method','matrix_exponential','normalize_by_mass',true, 'calculate_mass', true};
opts = OPTS(args);
pde = mirror1_collision_div(opts);
num_dims = numel(pde.dimensions);
% modify PDE
p.a.T_eV = 100;
p.a.v_beam = 0;
% run PDE
[err,fval,fval_realspace,nodes,err_realspace] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = err / norm(fval);
verifyLessThan(testCase,rel_err,1e-5);
end

function mirror1_collision_div_shiftBeam_test(testCase)
addpath(genpath(pwd));
disp('Testing the collision operator within mirror1 by shifting the beam');
% setup PDE
args = {'lev',5,'deg',4,'dt',1e-3,'quiet',true,'num_steps',10,'timestep_method','matrix_exponential','normalize_by_mass',true, 'calculate_mass', true};
opts = OPTS(args);
pde = mirror1_collision_div(opts);
num_dims = numel(pde.dimensions);
% modify PDE
p.a.E_eV = 4*pde.params.b.E_eV;
p.a.T_eV = pde.params.b.T_eV;
% run PDE
[err,fval,fval_realspace,nodes,err_realspace] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = err / norm(fval);
verifyLessThan(testCase,rel_err,1e-5);
end

function mirror1_collision_div_diffTemp_shiftBeam_test(testCase)
addpath(genpath(pwd));
disp('Testing the collision operator within mirror1 for different temperatures');
% setup PDE
args = {'lev',5,'deg',4,'dt',1e-3,'quiet',true,'num_steps',10,'timestep_method','matrix_exponential','normalize_by_mass',true, 'calculate_mass', true};
opts = OPTS(args);
pde = mirror1_collision_div(opts);
num_dims = numel(pde.dimensions);
% modify PDE
p.a.E_eV = 4*pde.params.b.E_eV;
p.a.T_eV = 4*pde.params.b.T_eV;
% run PDE
[err,fval,fval_realspace,nodes,err_realspace] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = err / norm(fval);
verifyLessThan(testCase,rel_err,1e-5);
end


function mirror2_advection_Efield_div_test(testCase)

addpath(genpath(pwd));
disp('Testing the 2D Electric field advection div equation');
% setup PDE  
args = {'lev',5,'deg',3,'dt',1e-6,'calculate_mass',true,'normalize_by_mass',true,'timestep_method','BE','quiet',true,'num_steps',5,'case',4};
opts = OPTS(args);
pde = mirror2_advection_Efield_div(opts);
num_dims = numel(pde.dimensions);

ic_v = @(x,p,t) x.^2;
ic1 = new_md_func(num_dims,{ic_v,[]});
pde.initial_conditions = {ic1};

src_v = @(v,p,t) v;
src_z = @(z,p,t) 2*pde.params.a.Z.*pde.params.e.*pde.params.E.*cos(z)./pde.params.a.m;
src = new_md_func(num_dims,{src_v,src_z});
pde.sources = {src};

    function ret = solution(v,p,t)
        ret =  v.^2;
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

function mirror1_space_div_test(testCase)
addpath(genpath(pwd));
disp('Testing the 1D spatial dimension');
% setup PDE
args = {'lev',5,'deg',5,'dt',1e-6,'quiet',true,'num_steps',5,'timestep_method','matrix_exponential','normalize_by_mass',true, 'calculate_mass', true};
opts = OPTS(args);
pde = mirror1_space_div(opts);
% num_dims = numel(pde.dimensions);
% ic1 = new_md_func(num_dims,{pde.sol_s});
% pde.initial_conditions = {ic1};
% run PDE
[err,fval,fval_realspace,nodes,err_realspace] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = err / norm(fval);
verifyLessThan(testCase,rel_err,1e-5);
end


function mirror1_space_diffSign_test(testCase)
addpath(genpath(pwd));
disp('Testing the 1D spatial dimension by changing the sign through the pitch angle');
% setup PDE
args = {'lev',5,'deg',5,'dt',1e-6,'quiet',true,'num_steps',5,'timestep_method','matrix_exponential','normalize_by_mass',true, 'calculate_mass', true};
opts = OPTS(args);
pde = mirror1_space_div(opts);

pde.params.pitch_test = 3*pi/4;
% num_dims = numel(pde.dimensions);
% ic1 = new_md_func(num_dims,{pde.sol_s});
% pde.initial_conditions = {ic1};
% run PDE
[err,fval,fval_realspace,nodes,err_realspace] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = err / norm(fval);
verifyLessThan(testCase,rel_err,1e-5);
end

function mirror1_space_diffSign2_test(testCase)
addpath(genpath(pwd));
disp('Testing the 1D spatial dimension by changing the sign through the magnetic field');
% setup PDE
args = {'lev',5,'deg',5,'dt',1e-6,'quiet',true,'num_steps',5,'timestep_method','matrix_exponential','normalize_by_mass',true, 'calculate_mass', true};
opts = OPTS(args);
pde = mirror1_space_div(opts);

pde.params.B_func = @(s) s;
pde.params.dB_ds = @(s) s.*0 + 1;
% num_dims = numel(pde.dimensions);
% ic1 = new_md_func(num_dims,{pde.sol_s});
% pde.initial_conditions = {ic1};
% run PDE
[err,fval,fval_realspace,nodes,err_realspace] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = err / norm(fval);
verifyLessThan(testCase,rel_err,1e-5);
end

function mirror2_pitch_space_div_test(testCase)

addpath(genpath(pwd));
disp('Testing the 2D pitch and spatial advection div equation');
% setup PDE  
args = {'lev',5,'deg',3,'dt',1e-6,'calculate_mass',true,'normalize_by_mass',false,'timestep_method','BE','quiet',true,'num_steps',3};
opts = OPTS(args);
pde = mirror2_space_pitch_div(opts);
num_dims = numel(pde.dimensions);

% modify PDE
ic_th = @(z,p,t,dat) cos(z);
ic_s = @(s,p,t,dat) s;
ic1 = new_md_func(num_dims,{ic_th,ic_s});
pde.initial_conditions = {ic1};

pde.params.B_func = @(s) s.^2;
pde.params.dB_ds = @(s) 2.*s;

src_th = @(x,p,t,d) (sin(x)).^2 - (cos(x)).^2;
src = new_md_func(num_dims,{src_th,[]});
pde.sources = {src};

    function ret = solution_th(z,p,t)
        ret =  cos(z);
        if isfield(p,'norm_fac')
            ret = p.norm_fac .* ret;
        end
    end


    function ret = solution_s(s,p,t)
        ret =  s;
        if isfield(p,'norm_fac')
            ret = p.norm_fac .* ret;
        end
    end

soln1 = new_md_func(num_dims,{@solution_th,@solution_s});
pde.solutions = {soln1};

% run PDE
[err,fval,fval_realspace,nodes,err_realspace,outputs] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = outputs.rel_err{end};
verifyLessThan(testCase,rel_err,1.5e-4);
end

% function mirror3_velocity_test(testCase)
% 
% addpath(genpath(pwd));
% disp('Testing the velocity dimension within mirror3');
% % setup PDE
% args = {'lev',3,'deg',3,'dt',10e-8,'calculate_mass',true,'normalize_by_mass',true,'timestep_method','matrix_exponential','quiet',true,'num_steps',1};
% opts = OPTS(args);
% pde = mirror3(opts);
% num_dims = numel(pde.dimensions);
% 
% % modify PDE
% ic_v = @(x,p,t) p.init_cond_v(x);
% ic1 = new_md_func(num_dims,{ic_v,[],[]});
% pde.initial_conditions = {ic1};
% 
% pde.params.boundary_cond_v = @(v,p,t) p.init_cond_v(v);
% pde.params.boundary_cond_z = @(z,p,t) z.*0 + 1;
% pde.params.boundary_cond_s = @(s,p,t) s.*0 + 1;
% 
% pde.dimensions{1,1}.min = 0;
% pde.dimensions{1,1}.max = 5e6;
% pde.dimensions{1,2}.min = 0;
% pde.dimensions{1,2}.max = pi/2;
% pde.terms = {pde.terms{1,1}, pde.terms{1,2}};
% 
%     function ret = solution(v,p,t)
%         ret =  p.soln_v(v,p,t);
%         if isfield(p,'norm_fac')
%             ret = p.norm_fac .* ret;
%         end
%     end
%     function ret = solution_z(z,p,~)
%         ret = 1/sin(z).^2;
%     end
% 
% soln1 = new_md_func(num_dims,{@solution,[],[]});
% pde.solutions = {soln1};
% 
% % run PDE
% [err,fval,fval_realspace,nodes,err_realspace,outputs] = asgard_run_pde(opts,pde);
% % assert on correctness
% rel_err = outputs.rel_err{end};
% verifyLessThan(testCase,rel_err,1.5e-4);
% end

% function mirror3_space_test(testCase)
% addpath(genpath(pwd));
% disp('Testing the spatial dimension within mirror3');
% % setup PDE
% args = {'lev',3,'deg',3,'dt',5e-15, 'normalize_by_mass', true, 'timestep_method', 'matrix_exponential','quiet',true,'num_steps',1};
% opts = OPTS(args);
% pde = mirror3(opts);
% num_dims = numel(pde.dimensions);
% 
% % modify PDE
% ic_s = @(x,p,t) exp(x);
% ic1 = new_md_func(num_dims,{[],[],ic_s});
% pde.initial_conditions = {ic1};
% 
% pde.params.boundary_cond_v = @(v,p,t) v.*0 + 1;
% pde.params.boundary_cond_z = @(z,p,t) z.*0 + 1;
% pde.terms = {pde.terms{1,4}, pde.terms{1,5}};
% 
%     function ret = solution(s,p,t)
%         ret =  p.soln_s(s,p,t);
%         if isfield(p,'norm_fac')
%             ret = p.norm_fac .* ret;
%         end
%     end
% 
% soln1 = new_md_func(num_dims,{[],[],@solution});
% pde.solutions = {soln1};
% 
% % run PDE
% [err,fval,fval_realspace,nodes,err_realspace,outputs] = asgard_run_pde(opts,pde);
% % assert on correctness
% rel_err = outputs.rel_err{end};
% verifyLessThan(testCase,rel_err,1e-4);
% end
% 
% function mirror3_manufacture_test(testCase)
% addpath(genpath(pwd));
% disp('Testing a manufactured solution within mirror3');
% % setup PDE
% args = {'lev',3,'deg',3,'dt',1e-12,'calculate_mass',true,'quiet',true,'num_steps',1,'normalize_by_mass',true,'timestep_method','BE'};
% opts = OPTS(args);
% pde = mirror3(opts);
% num_dims = numel(pde.dimensions);
% 
% % modify PDE
% pde.params.nu_D = @(v,a,b) v;
% pde.params.nu_s = @(v,a,b) v;
% pde.params.nu_par = @(v,a,b) v;
% pde.params.boundary_cond_v = @(v,p,t) exp(-v.*t);
% pde.params.boundary_cond_s = @(s,p,t) s.*0 + 1;
% pde.terms = {pde.terms{1,1},pde.terms{1,2},pde.terms{1,3},pde.terms{1,4}};
% 
% %source function for manufactured solution
% src_v = @(v,p,t) exp(-v.*t).*((p.a.m./(p.a.m + p.b.m)).*(-4*v + v.^2.*t) + 5*v.^2.*t./2 - v.^3.*t.^2./2);
% src_z = @(z,p,t) cos(z);
% 
% s_v = @(v,p,t) exp(-v.*t);
% s_z = @(z,p,t) cos(z);
% soln1 = new_md_func(num_dims,{s_v,s_z,[]});
% src = new_md_func(num_dims,{src_v,src_z,[]});
% pde.solutions = {soln1};
% pde.sources = {src};
% 
% ic_z = @(z,p,t) cos(z);
% ic1 = new_md_func(num_dims,{[],ic_z,[]});
% pde.initial_conditions = {ic1};
% 
% % run PDE
% [err,fval,fval_realspace,nodes,err_realspace,outputs] = asgard_run_pde(opts,pde);
% % assert on correctness
% rel_err = outputs.rel_err{end};
% verifyLessThan(testCase,rel_err,1e-4);
% end