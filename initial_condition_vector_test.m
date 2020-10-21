function tests = initial_condition_vector_test()
fh = localfunctions;
tests = functiontests(fh);
end

function sane_1D_function_test(testCase)
addpath(genpath(pwd));

opts = OPTS();
opts.lev = 4;
opts.deg = 4;
opts.max_lev_coeffs = true;
pde = advection1(opts);
num_dims = numel(pde.dimensions);

    function ans = my_func(x)
        ans = cos(x);
%         ans = 1./(x.^2);
    end

pde.dimensions{1}.min = 0.1;
pde.dimensions{1}.max = 2*pi;
pde.dimensions{1}.init_cond_fn = @(x,p,t) my_func(x);

[elements, elements_idx]    = hash_table_sparse_nD (pde.get_lev_vec, opts.max_lev, opts.grid_type);
hash_table.elements         = elements;
hash_table.elements_idx     = elements_idx;

for d=1:num_dims
[Meval{d},nodes{d}] = matrix_plot_D(pde,opts,pde.dimensions{d});
end

t = 0;
fval = initial_condition_vector(pde,opts,hash_table,t);
fval_realspace = wavelet_to_realspace(pde,opts,Meval,fval,hash_table); 

do_plot = 0;
if do_plot
plot(nodes{1},fval_realspace,'-o')
hold on
plot(nodes{1},my_func(nodes{1}),'-o')
hold off
end
err = norm(fval_realspace - my_func(nodes{1}'));
verifyLessThan(testCase, err, 4.4e-4);
end
