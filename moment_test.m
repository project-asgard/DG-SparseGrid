function tests = moment_test()
tests = functiontests(localfunctions);
end

function test_moment_2D_spherical(testCase)

args = {'lev',5,'deg',4,'case',3};
opts = OPTS(args);
pde = diffusion2(opts);

num_dims = numel(pde.dimensions);
[elements, elements_idx]    = hash_table_sparse_nD (pde.get_lev_vec, opts.max_lev, opts.grid_type);
hash_table.elements         = elements;
hash_table.elements_idx     = elements_idx;
t = 0;

%% Verify moment f(x,y) = 1 (should be auto added)

ic_x = @(x,p,t) cos(pi*x);
ic_y = @(y,p,t) cos(pi*y);
ic_t = @(t,p) exp(-2*pi^2*t);
pde.initial_conditions = {new_md_func(num_dims,{ic_x,ic_y,ic_t})};

pde = compute_dimension_mass_mat(opts,pde);
pde = calculate_moment_data(pde,opts);

fval = initial_condition_vector(pde,opts,hash_table,t);
[pde, moment_val, ~] = calculate_mass(pde,opts,fval,hash_table,t);

%moment_val should be 0
diff1 = abs(moment_val-0);
verifyLessThan(testCase, diff1, 1e-10);

%% Verify moment f(x,y) = ic

pde.moments{1} = MOMENT({new_md_func(num_dims,{ic_x,ic_y,ic_t})});

pde = compute_dimension_mass_mat(opts,pde);
pde = calculate_moment_data(pde,opts);

fval = initial_condition_vector(pde,opts,hash_table,t);
[pde, moment_val, ~] = calculate_mass(pde,opts,fval,hash_table,t);

%moment_val should be 0
diff1 = abs(moment_val-0.25);
verifyLessThan(testCase, diff1, 1e-13);

%% Verify moment f(x,y) = y

ic_x = @(x,p,t) 1-x.^2;
ic_y = @(y,p,t) y;
ic_t = @(t,p) 0*t+1;
pde.initial_conditions = {new_md_func(num_dims,{ic_x,ic_y,ic_t})};

mom_x = @(x,p,t) 0*x+1;
mom_y = @(y,p,t) y;
mom_t = @(t,p) 0*t+1;
mom_func = new_md_func(num_dims,{mom_x,mom_y,mom_t});

pde.moments{1} = MOMENT({mom_func});

pde = compute_dimension_mass_mat(opts,pde);
pde = calculate_moment_data(pde,opts);

fval = initial_condition_vector(pde,opts,hash_table,t);
[pde, moment_val, ~] = calculate_mass(pde,opts,fval,hash_table,t);

%moment_val should be 2/9
diff1 = abs(moment_val-2/9);
verifyLessThan(testCase, diff1, 1e-13);


%% Do non-cartesian coordinate example

pde = diffusion2_spherical(opts);

num_dims = numel(pde.dimensions);
[elements, elements_idx]    = hash_table_sparse_nD (pde.get_lev_vec, opts.max_lev, opts.grid_type);
hash_table.elements         = elements;
hash_table.elements_idx     = elements_idx;
t = 0;

mom_x = @(r,p,t) 0*r+1;
mom_y = @(th,p,t) cos(th);
mom_t = @(t,p) 0*t+1;
mom_func = new_md_func(num_dims,{mom_x,mom_y,mom_t});
pde.moments{1} = MOMENT({mom_func});

pde = compute_dimension_mass_mat(opts,pde);
pde = calculate_moment_data(pde,opts);

fval = initial_condition_vector(pde,opts,hash_table,t);
[pde, moment_val, ~] = calculate_mass(pde,opts,fval,hash_table,t);

%moment_val should be 0
diff1 = abs(moment_val-8/3);
verifyLessThan(testCase, diff1, 1e-13);


end