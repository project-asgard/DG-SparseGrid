function tests = mirror_test()

fh = localfunctions;
tests = functiontests(fh);

end

function mirror1_collision_sameTemp_test(testCase)
addpath(genpath(pwd));
disp('Testing the collision operator within mirror1');
% setup PDE
args = {'lev',5,'deg',4,'dt',1e-3,'quiet',true,'num_steps',2,'timestep_method','matrix_exponential','normalize_by_mass',true, 'calculate_mass', true};
opts = OPTS(args);
pde = mirror1_collision(opts);
num_dims = numel(pde.dimensions);
% modify PDE
ic_v = @(v,p,t) p.soln_v(v);
ic1 = new_md_func(num_dims,{ic_v});
pde.initial_conditions = {ic1};
% run PDE
[err,fval,fval_realspace,nodes,err_realspace] = asgard_run_pde(opts,pde);
% assert on correctness
rel_err = err / norm(fval);
verifyLessThan(testCase,rel_err,1e-5);
end

function mirror1_collision_diffTemp_test(testCase)
addpath(genpath(pwd));
disp('Testing the collision operator within mirror1 by changing the temperature of the beam');
% setup PDE
args = {'lev',5,'deg',4,'dt',1e-3,'quiet',true,'num_steps',10,'timestep_method','matrix_exponential','normalize_by_mass',true, 'calculate_mass', true};
opts = OPTS(args);
pde = mirror1_collision(opts);
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

function mirror1_collision_shiftBeam_test(testCase)
addpath(genpath(pwd));
disp('Testing the collision operator within mirror1 by shifting the beam');
% setup PDE
args = {'lev',5,'deg',4,'dt',1e-3,'quiet',false,'num_steps',10,'timestep_method','matrix_exponential','normalize_by_mass',true, 'calculate_mass', true};
opts = OPTS(args);
pde = mirror1_collision(opts);
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

function mirror1_collision_diffTemp_shiftBeam_test(testCase)
addpath(genpath(pwd));
disp('Testing the collision operator within mirror1');
% setup PDE
args = {'lev',5,'deg',4,'dt',1e-3,'quiet',false,'num_steps',10,'timestep_method','matrix_exponential','normalize_by_mass',true, 'calculate_mass', true};
opts = OPTS(args);
pde = mirror1_collision(opts);
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

function mirror1_space_div_test(testCase)

addpath(genpath(pwd));
disp('Testing the spatial dimension within mirror1');
% setup PDE
args = {'lev',3,'deg',3,'dt',5e-7,'quiet',true,'num_steps',20,'timestep_method','matrix_exponential','normalize_by_mass',false};
opts = OPTS(args);
pde = mirror1_space_div(opts);

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

function mirror2_advection_Efield_div_test(testCase)

addpath(genpath(pwd));
disp('Testing the 2D Electric field advection div equation');
% setup PDE  
args = {'lev',5,'deg',3,'dt',1e-5,'calculate_mass',true,'normalize_by_mass',true,'timestep_method','BE','quiet',false,'num_steps',5};
opts = OPTS(args);
pde = mirror2_advection_Efield_div(opts);
num_dims = numel(pde.dimensions);

% modify PDE
ic_v = @(x,p,t) x.^2;
ic1 = new_md_func(num_dims,{ic_v,[]});
pde.initial_conditions = {ic1};

n_cgs = 5e13; %equilibrium density in cm.^-3
m_e_cgs = 9.109*10^-28; %electron mass in g
m_D_cgs = 3.3443*10^-24; %Deuterium mass in g
m_He_cgs = 6.7*10^-24; %helium 4 mass in g 
m_B_cgs = 1.82*10^-23; %Boron 11 mass in g
temp_cgs = 1.6022e-10; %temperature in erg
params_cgs.a.n = n_cgs;
params_cgs.b.n = n_cgs;
params_cgs.b2.n = n_cgs;
params_cgs.a.m = m_e_cgs; %beam is electrons
params_cgs.b.m = m_e_cgs; %background is electrons
params_cgs.b2.m = m_B_cgs;
params_cgs.a.Z = -1;
params_cgs.b.Z = -1;
params_cgs.b2.Z = 1;
params_cgs.e = 4.803*10^-10; %charge in Fr
params_cgs.E = 2.6e-5; %E field in statvolt/cm
pde.params.a.n = 10^6*params_cgs.a.n;%converting to m^-3
pde.params.b.n = 10^6*params_cgs.b.n;
pde.params.b2.n = 10^6*params_cgs.b2.n;

pde.params.a.m = 0.001*params_cgs.a.m; %converting to kg
pde.params.b.m = 0.001*params_cgs.b.m; 
pde.params.b2.m = 0.001*params_cgs.b2.m;
pde.params.a.Z = params_cgs.a.Z;
pde.params.b.Z = params_cgs.b.Z;
pde.params.b2.Z = params_cgs.b2.Z;
pde.params.ln_delt = 27.0857;
pde.params.a.T_eV = 1.33e4;%2/3*params.a.E_eV;
pde.params.b.T_eV = pde.params.a.T_eV;
pde.params.b2.T_eV = pde.params.a.T_eV;
pde.params.a.vth = pde.params.v_th(pde.params.a.T_eV,pde.params.a.m);
pde.params.b.vth = pde.params.v_th(pde.params.b.T_eV,pde.params.b.m);
pde.params.b2.vth = pde.params.v_th(pde.params.b2.T_eV,pde.params.b2.m);
E_dreicer_si = pde.params.a.n.*pde.params.e^3*pde.params.ln_delt/(2*pi*pde.params.eps0^2*pde.params.a.m ... 
            *pde.params.a.vth^2);
pde.params.E = 10^-2*E_dreicer_si; 
src_v = @(v,p,t) v;
src_z = @(z,p,t) 2*pde.params.a.Z.*pde.params.e.*pde.params.E.*cos(z)./pde.params.a.m;
src = new_md_func(num_dims,{src_v,src_z});
pde.sources = {src};

pde.params.boundary_cond_v = @(v,p,t) v.^2;
pde.params.boundary_cond_z = @(z,p,t) z.*0 + 1;
pde.params.boundary_cond_s = @(s,p,t) s.*0 + 1;

pde.dimensions{1,1}.min = 0;
pde.dimensions{1,1}.max = 1e9;
pde.dimensions{1,2}.min = 0;
pde.dimensions{1,2}.max = pi;
pde.terms = {pde.terms{1,1}, pde.terms{1,2},pde.terms{1,3}};

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
    function ret = solution_z(z,p,~)
        ret = 1/sin(z).^2;
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
