function tests = moment_integral_test()
tests = functiontests(localfunctions);
end

function test_moment_2D_spherical(testCase)

args = {'lev',4,'deg',4,'case',3};
opts = OPTS(args);
pde = mirror_velocity2(opts);

num_dims = numel(pde.dimensions);

    function ans = my_func_v(v)
        v_th = 0.2;
        ans = exp(-v.^2./v_th^2);
    end
    function ans = my_func_z(z)
        ans = z.*+1;
    end

pde.dimensions{1}.min = 0;
pde.dimensions{1}.max = 1;
pde.dimensions{1}.moment_dV = @(v,p,t) v.^2;

pde.dimensions{2}.min = 0;
pde.dimensions{2}.max = pi;
pde.dimensions{2}.moment_dV = @(z,p,t) sin(z);

pde = compute_dimension_mass_mat(opts,pde);

% overwrite the initial conditions for this test

ic_v = @(v,p,t) my_func_v(v);
ic_z = @(z,p,t) my_func_z(z);
ic1 = new_md_func(num_dims,{ic_v,ic_z});
pde.initial_conditions = {ic1};

[elements, elements_idx]    = hash_table_sparse_nD (pde.get_lev_vec, opts.max_lev, opts.grid_type);
hash_table.elements         = elements;
hash_table.elements_idx     = elements_idx;

for d=1:num_dims
[Meval{d},nodes{d}] = matrix_plot_D(pde,opts,pde.dimensions{d});
end
coord = get_realspace_coords(pde,nodes);

t = 0;
fval = initial_condition_vector(pde,opts,hash_table,t);
fval_realspace = wavelet_to_realspace(pde,opts,Meval,fval,hash_table); 

do_plot = 0;
if do_plot
    hold off;
    contour(nodes{1},nodes{2},reshape(fval_realspace,numel(nodes{1}),numel(nodes{2})));
end

% check mass moment

    function ans = mass_func_v(v,p,t)
        ans = v.*0+1;
    end
    function ans = mass_func_z(z,p,t)
        ans = z.*0+1;
    end

moment_func_nD = {@mass_func_v,@mass_func_z};
mass = moment_integral(opts.lev, opts.deg, coord, fval_realspace, moment_func_nD, pde.dimensions);
mass_analytic = 0.0111367;
diff1 = abs(mass-mass_analytic);
verifyLessThan(testCase, diff1, 1e-6);

% check energy moment
    function ans = energy_func_v(v,p,t)
        ans = v.^2;
    end

moment_func_nD = {@energy_func_v,@mass_func_z};
energy = moment_integral(opts.lev, opts.deg, coord, fval_realspace, moment_func_nD, pde.dimensions);
energy_analytic = 0.000668199;
diff2 = abs(energy-energy_analytic);
verifyLessThan(testCase, diff2, 1e-6);

% check some 2D moment
    function ans = random_func_v(v,p,t)
        ans = v.^2;
    end
    function ans = random_func_z(z,p,t)
        ans = z.^2;
    end

moment_func_nD = {@random_func_v,@random_func_z};
moment = moment_integral(opts.lev, opts.deg, coord, fval_realspace, moment_func_nD, pde.dimensions);
moment_analytic = 0.00258567;
diff3 = abs(moment-moment_analytic);
verifyLessThan(testCase, diff3, 1e-6);

end
